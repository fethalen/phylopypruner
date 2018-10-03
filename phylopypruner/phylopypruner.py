#!/usr/bin/env python
#pylint: disable=too-many-branches
#pylint: disable=too-many-locals

# Author: Felix Thalen
# Date: July 24, 2018

r"""
  PhyloPyPruner: modular and contamination aware tree-based orthology inference

  Refer to the wiki (https://gitlab.com/fethalen/phylopypruner/wikis/home) for a
  tutorial and explanations of the implemented algorithms.

example:
  ./phylopypruner.py --msa 16s.fa --tree 16s.tre --min-taxa 4 --min-len 200
    --min-support 0.7 --trim-lb 5 --outgroup Drosophila --root midpoint
    --mask pdist --prune MI

"""

from __future__ import print_function
import argparse
import os
import sys
import datetime
import shutil
from textwrap import wrap
import fasta
import newick
import filtering
import decontamination
import mask_monophylies
import root
from prune_paralogs import prune_paralogs
from summary import Summary
from summary import mk_sum_out_title
from log import Log
from settings import Settings
# ensures that input is working across both Python 2 and 3
if hasattr(__builtins__, "raw_input"): input = raw_input

VERSION = 0.1
FASTA_EXTENSIONS = {".fa", ".fas", ".fasta", ".fna", ".faa", ".fsa", ".ffn",
                    ".frn"}
NW_EXTENSIONS = {".newick", ".nw", ".tre", ".tree"}
HEADER = "id;sequences;seq_len_avg;shortest_seq;longest_seq;pct_missing_data;\
alignment_len\n"
TIMESTAMP = datetime.datetime.now().strftime("%Y-%m-%d")
ORTHO_STATS_PATH = "/{}_ppp_ortho_stats.csv".format(TIMESTAMP)
LOG_PATH = "/{}_ppp_run.log".format(TIMESTAMP)
ORTHOLOGS_PATH = "/{}_orthologs".format(TIMESTAMP)

def _warning(message):
    """
    Returns the provided message with the text 'warning: ' in bold red
    prepended.
    """
    return "\n".join(wrap("{}warning{}: {}".format("\033[91m\033[1m",
                                                   "\033[0m", message), 80))

NO_FILES = _warning("""Please provide either a multiple sequence alignment
(MSA) and a Newick tree, or a path to a directory containing multiple MSAs
and Newick trees. Run PhyloPyPruner using the '-h' or '--help' flag for
additional instructions.""")

def _yes_or_no(question):
    """
    Takes a question as an input and prompts the user for a yes or no. Returns
    True if the answer is yes and False if the answer is no.
    """
    # make input work the same way in both Python 2 and 3
    answer = input(question + " (y/n): ".lower().rstrip())
    while not (answer == "y" or answer == "n"):
        print("type 'y' for yes and 'n' for no")
        answer = input(question + " (y/n): ".lower().rstrip())
    if answer[0] == "y":
        return True
    elif answer[0] == "n":
        return False

def _validate_input(msa, tree, tree_path):
    "Test to see if MSA and tree entries matches."
    descriptions = list(msa.iter_descriptions())
    names = list(tree.iter_names())

    if set(descriptions).intersection(names) < set(descriptions):
        raise AssertionError("MSA {} does not match tree \
                {}".format(msa.filename, tree_path))

def _no_files(args):
    "Returns true if the required files are not provided."
    if not args.msa and not args.tree and not args.dir:
        return True
    elif args.msa and not args.tree or not args.msa and args.tree:
        return True
    else:
        return False

def _validate_arguments(args):
    "Performs a series of checks to validate the input before execution."
    valid_arguments = True

    if _no_files(args):
        print(NO_FILES)
        valid_arguments = False

    if not args.outgroup and args.prune == "MO" or \
        not args.outgroup and args.prune == "RT":
        print(_warning("trying to run {}, but no outgroup has been\
 specified".format(args.prune)))
        valid_arguments = False

    if args.min_support:
        if args.min_support < 0 or args.min_support > 100:
            print(_warning("minimum support value ('--min-support') has to be\
 either in percentage (1-100) or in decimal format between 0.0 - 1.0"))
            valid_arguments = False
        elif args.min_support > 1:
            # convert from percentage to floating point
            args.min_support = args.min_support / 100

    if args.min_len and args.min_len < 1:
        print(_warning("minimum sequence length ('--min-len') has to be a\
 positive integer (1, 2, 3, 4, ...)"))
        valid_arguments = False

    if args.min_taxa and args.min_taxa < 1:
        print(_warning("minimum number of taxa ('--min-taxa') has to be a\
 positive integer (1, 2, 3, 4, ...)"))
        valid_arguments = False

    if args.trim_lb and args.trim_lb <= 0:
        print(_warning("the factor for removing long branches ('--trim-lb')\
has to be a positive number"))
        valid_arguments = False

    if not valid_arguments:
        exit()

def _run(settings, msa, tree):
    """
    Takes a dictionary that specifies a set of settings as an input. The
    dictionary should contain the following keys: msa, tree, min_taxa, min_len,
    min_support, trim_lb, outgroup, root, mask, and prune. Runs PhyloPyPruner
    once using these settings.
    """
    # test to see if the entries in the MSA and tree matches
    _validate_input(msa, tree, settings.nw_file)

    log = Log(VERSION, msa, tree, settings)
    log.homology_tree = tree.view()

    # remove short sequences
    if settings.min_len:
        log.trimmed_seqs = filtering.trim_short_seqs(msa, tree, settings.min_len)

    # exit if number of OTUs < threshold
    if settings.min_taxa:
        if filtering.too_few_otus(tree, settings.min_taxa):
            return log

    # trim long branches
    if settings.trim_lb:
        log.lbs_removed = list(filtering.prune_long_branches(tree, settings.trim_lb))

    # exit if number of OTUs < threshold
    if settings.min_taxa:
        if filtering.too_few_otus(tree, settings.min_taxa):
            return log

    # collapse weakly supported nodes into polytomies
    if settings.min_support:
        log.pruned_sequences = filtering.collapse_nodes(tree, settings.min_support)

    # mask monophyletic groups
    if settings.mask:
        if settings.mask == "pdist":
            tree, masked_seqs = mask_monophylies.pairwise_distance(tree)
        elif settings.mask == "longest":
            tree, masked_seqs = mask_monophylies.longest_isoform(msa, tree)
        log.monophylies_masked = masked_seqs

    # exit if number of OTUs < threshold
    if settings.min_taxa:
        if filtering.too_few_otus(tree, settings.min_taxa):
            return log

    # root by outgroup
    if settings.outgroup:
        if not settings.prune == "MO":
            tree = root.outgroup(tree, settings.outgroup)

    # exclude taxa within the list settings.exclude
    if settings.exclude:
        tree = filtering.exclude(tree, settings.exclude)

    # exit if number of OTUs < threshold
    if settings.min_taxa:
        if filtering.too_few_otus(tree, settings.min_taxa):
            return log

    # mask monophyletic groups
    if settings.mask:
        if settings.mask == "pdist":
            tree, masked_seqs = mask_monophylies.pairwise_distance(tree)
        elif settings.mask == "longest":
            tree, masked_seqs = mask_monophylies.longest_isoform(msa, tree)
        log.monophylies_masked.update(masked_seqs)

    log.masked_tree = tree
    log.masked_tree_str = tree.view()

    # exit if number of OTUs < threshold
    if settings.min_taxa:
        if filtering.too_few_otus(tree, settings.min_taxa):
            return log

    # get a list of paralogs
    log.paralogs = tree.paralogs()

    # prune paralogs
    log.orthologs = prune_paralogs(settings.prune, tree,
                                   settings.min_taxa, settings.outgroup)

    return log

def _get_orthologs(settings, directory="", dir_out=None, verbose=False):
    fasta_path = "{}{}".format(directory, settings.fasta_file)
    nw_path = "{}{}".format(directory, settings.nw_file)
    msa = fasta.read(fasta_path)
    nw_file = newick.read(nw_path)
    log = _run(settings, msa, nw_file)

    log.get_msas_out(dir_out)
    log.report(verbose, dir_out)
    return log

def parse_args():
    """
    Parse the arguments provided by the user.
    """
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("-v", "--version",
                        action="version",
                        version=str(VERSION),
                        help="display the version number and exit")
    parser.add_argument("--msa",
                        metavar="<MSA>",
                        type=str,
                        help="multiple sequence alignment (MSA) in FASTA\
                              format")
    parser.add_argument("--tree",
                        metavar="<tree>",
                        type=str,
                        help="Newick tree file")
    parser.add_argument("--dir",
                        metavar="<directory>",
                        type=str,
                        default=None,
                        help="directory containing multiple MSAs and trees\
                              for analysis")
    parser.add_argument("--output",
                        metavar="<directory>",
                        type=str,
                        default=None,
                        help="path to output directory")
    parser.add_argument("--min-taxa",
                        metavar="<threshold>",
                        type=int,
                        default=None,
                        help="minimum number of OTUs allowed in output")
    parser.add_argument("--min-len",
                        metavar="<threshold>",
                        default=None,
                        type=int,
                        help="remove sequences shorter than <threshold>")
    parser.add_argument("--min-support",
                        metavar="<threshold>",
                        default=None,
                        type=float,
                        help="collapse nodes with a support value below\
                              <threshold> into polytomies; in percentage\
                              or decimal format")
    parser.add_argument("--trim-lb",
                        default=None,
                        metavar="<factor>",
                        type=float,
                        help="remove sequences with a branch length that is\
                              <factor> times longer than the standard\
                              deviation of all branches")
    parser.add_argument("--outgroup",
                        nargs='+',
                        metavar="<OTU>",
                        default=None,
                        type=str,
                        help="one or more outgroup OTUs; root using outgroup\
                              rooting, if the outgroup forms a non-repetetive,\
                              monophyletic group")
    parser.add_argument("--root",
                        default=None,
                        type=str,
                        choices=["midpoint", "molecular_clock"],
                        help="reroot tree using this rooting method if an\
                              outgroup hasn't been provided or if outgroup\
                              rooting fails")
    parser.add_argument("--mask",
                        default="pdist",
                        type=str,
                        choices=["longest", "pdist"],
                        help="mask monophylies using this method")
    parser.add_argument("--prune",
                        default="LS",
                        type=str,
                        choices=["LS", "MI", "MO", "RT", "1to1"],
                        help="prune paralogs using this method")
    parser.add_argument("--exclude",
                        nargs='+',
                        metavar="<OTU>",
                        default=None,
                        type=str,
                        help="a list of OTUs to exclude in this run")
    parser.add_argument("--include",
                        nargs='+',
                        metavar="<OTU>",
                        default=None,
                        type=str,
                        help="always include these OTUs, even if they were \
                        flagged as problematic by '--trim-freq-paralogs' or \
                        '--trim-divergent'")
    parser.add_argument("--jackknife",
                        default=False,
                        action="store_true",
                        help="leave out OTUs one by one and save statistics \
                        for each OTU removed in the summary; note that no \
                        data is produced, use the --exclude flag to exclude \
                        OTUs in a subsequent analysis")
    parser.add_argument("--trim-freq-paralogs",
                        default=None,
                        metavar="<factor>",
                        type=float,
                        help="remove OTUs with a paralogy frequency that is \
                              than <factor> times more frequent than the \
                              standard deviation of all OTUs together")
    parser.add_argument("--trim-divergent",
                        default=None,
                        metavar="<factor>",
                        type=float,
                        help="remove OTUs where the average pairwise distance \
                              of the sequences for that OTU is larger than \
                              <factor> times the average pairwise distance of \
                              all sequences")
    parser.add_argument("--wrap",
                        metavar="<max column>",
                        default=None,
                        type=int,
                        help="wrap output sequences at column <max column>")
    parser.add_argument("--verbose",
                        default=False,
                        action="store_true",
                        help="show a more detailed report")
    return parser.parse_args(args=None if sys.argv[1:] else ['--help'])

def main():
    "Parse args, run filter and infer orthologs."
    args = parse_args()
    _validate_arguments(args)
    settings = Settings(args)
    summary = Summary()

    dir_out = None
    if args.output:
        # create output directory if it does not currently exist
        dir_out = args.output.rstrip("/") # get rid of trailing slash in path
    elif args.dir:
        dir_out = str(args.dir).rstrip("/")
    else:
        dir_out = os.path.dirname(str(args.msa).rstrip("/"))
    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)

    if os.path.isfile(dir_out + ORTHO_STATS_PATH):
        question = _warning("files from a previous run exists in the output \
directory, overwrite?")
        if not _yes_or_no(question):
            exit()

    if os.path.isfile(dir_out + ORTHO_STATS_PATH):
        os.remove(dir_out + ORTHO_STATS_PATH)
    with open(dir_out + ORTHO_STATS_PATH, "w") as stats_file:
        stats_file.write(HEADER)

    if os.path.isfile(dir_out + LOG_PATH):
        os.remove(dir_out + LOG_PATH)

    if os.path.isdir(dir_out + ORTHOLOGS_PATH):
        shutil.rmtree(dir_out + ORTHOLOGS_PATH)

    print("PhyloPyPruner version {}".format(VERSION))
    print(datetime.datetime.now().strftime("%A, %d. %B %Y %I:%M%p"))

    if args.msa and args.tree:
        # run for a single pair of files
        summary.logs.append(_get_orthologs(settings, "", dir_out,
                                           args.verbose))
    elif args.dir:
        # run for multiple files in directory
        if not args.dir[-1] == "/":
            dir_in = args.dir + "/"
        else:
            dir_in = args.dir
        if not os.path.isdir(dir_in):
            print(_warning("input directory {} does not exist".format(dir_in)))
            exit()
        # corresponding files; filename is key, values are tuple pairs where
        # MSA comes first and tree second
        corr_files = dict()
        for file in os.listdir(args.dir):
            filename, extension = os.path.splitext(file)
            extension = extension.lower()
            if extension in FASTA_EXTENSIONS or extension in NW_EXTENSIONS:
                if filename in corr_files.keys():
                    if extension in NW_EXTENSIONS:
                        corr_files[filename] = corr_files[filename] + (file,)
                    else:
                        corr_files[filename] = (file,) + corr_files[filename]
                else:
                    corr_files[filename] = (file,)

        total = len(corr_files)
        if total < 1:
            # no file pairs found in the provided directory
            print(_warning("no file pairs were found in the provided\
 directory"))
            exit()

        for index, pair in enumerate(corr_files, 1):
            try:
                settings.fasta_file, settings.nw_file = corr_files[pair]
            except:
                # corresponding file not found
                continue
            sys.stdout.flush()
            print("processing MSA: {}; processing tree: {} ({}/{} file \
pairs)".format(settings.fasta_file, settings.nw_file, index, total),
                  end="\r")
            summary.logs.append(_get_orthologs(settings, dir_in, dir_out,
                                               args.verbose))
        print("")
        mk_sum_out_title(dir_out)
        homolog_report = summary.homolog_report(dir_out)
        ortholog_report = summary.report("orthologs", dir_out)
        paralog_freq = summary.paralogy_frequency(dir_out)

    otus_to_exclude = []

    if args.trim_divergent:
        otus_to_exclude += decontamination.trim_divergent_otus(
            summary, args.trim_divergent)

    if args.trim_freq_paralogs:
        otus_to_exclude += decontamination.trim_freq_paralogs(
            summary, args.trim_freq_paralogs, paralog_freq)

    if otus_to_exclude:
        if args.include:
            otus_to_exclude = [otu for otu in otus_to_exclude if not otu in
                               args.include]
        summary, ortholog_report = decontamination.prune_by_exclusion(
            summary, otus_to_exclude, dir_out)

    if args.jackknife:
        decontamination.jackknife(summary, dir_out)

    print("{}\n{}".format(homolog_report, ortholog_report))

    summary.write_msas(args.wrap)

if __name__ == "__main__":
    main()
