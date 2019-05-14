#!/usr/bin/env python

# Author: Felix Thalen
# Date: July 24, 2018

"""
PhyloPyPruner is a Python package for tree-based orthology inference that is
used to refine the output of a graph-based approach by removing sequences
related via gene duplication. See gitlab.com/fethalen/phylopypruner for details.
"""

from __future__ import print_function
from __future__ import absolute_import
import argparse
import os
import shutil
import sys
import time
import datetime
from functools import partial
from multiprocessing import Pool
from multiprocessing import cpu_count
import pkg_resources
from phylopypruner import fasta
from phylopypruner import newick
from phylopypruner import filtering
from phylopypruner import decontamination
from phylopypruner import mask_monophylies
from phylopypruner import root
from phylopypruner import taxonomic_groups
from phylopypruner.prune_paralogs import prune_paralogs
from phylopypruner.summary import Summary
from phylopypruner.summary import mk_sum_out_title
from phylopypruner.log import Log
from phylopypruner.settings import Settings
from phylopypruner.supermatrix import Supermatrix
# ensures that input is working across Python 2.7 and 3+
if hasattr(__builtins__, "raw_input"): input = raw_input

VERSION = pkg_resources.require("phylopypruner")[0].version
FASTA_EXTENSIONS = {".fa", ".fas", ".fasta", ".fna", ".faa", ".fsa", ".ffn",
                    ".frn"}
NW_EXTENSIONS = {".newick", ".nw", ".tre", ".tree", ".out", ".treefile"}
ORTHO_STATS_PATH = "/output_alignment_stats.csv"
HOMOLOG_STATS_PATH = "/input_alignment_stats.csv"
ORTHOLOG_STATS_HEADER = "filename;otus;sequences;meanSeqLen;shortestSeq;longestSeq;\
pctMissingData;alignmentLen\n"
HOMOLOG_STATS_HEADER = "filename;otus;sequences;meanSeqLen;shortestSeq;longestSeq;\
pctMissingData;alignmentLen;shortSequencesRemoved;longBranchesRemoved;\
monophyliesMasked;nodesCollapsed;divergentOtusRemoved\n"
ORTHO_STATS_PATH = "/output_alignment_stats.csv"
LOG_PATH = "/phylopypruner.log"
ORTHOLOGS_PATH = "/output_alignments"
TIMESTAMP = datetime.datetime.now().strftime("%A, %d. %B %Y %I:%M%p")
ABOUT = "PhyloPyPruner version {}\n{}\n{}".format(
    VERSION, TIMESTAMP, "-" * len(TIMESTAMP))
NO_FILES = """Please provide either a multiple sequence alignment (MSA) and a
Newick tree, or a path to a directory containing multiple MSAs and Newick
trees."""

def _warning(message):
    """Print the provided message with the text 'warning: ' in front of it.

    Parameters
    ----------
    message : str
        Print this string

    Returns
    -------
    None
    """
    return "{}warning: {}{}".format(
        "\033[35m", "\033[0m", message)

def _error(message):
    """Print the provided message with the text 'error: ' in front of it and
    exit.

    Parameters
    ----------
    message : str
        Print this message before exiting.

    Returns
    -------
    None
    """
    print("{}error: {}{}".format("\033[31m", "\033[0m", message))
    exit()

def _yes_or_no(question):
    """
    Takes a question as an input and prompts the user for a yes or no. Returns
    True if the answer is yes and False if the answer is no.

    Parameters
    ----------
    question : str
        Print this message before (y/n).

    Returns
    -------
    True or False
        True if answer is 'y', False if answer is 'n'.
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
        print("example tree names:",names[:2], file=sys.stderr)
        print("example sequences:",descriptions[:2], file=sys.stderr)
        _error("MSA names don't match tree \n   {}\n   {}".format(msa.filename, tree_path))

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

    if _no_files(args):
        _error(NO_FILES)

    if not args.outgroup and args.prune == "MO" or \
        not args.outgroup and args.prune == "RT":
        _error("pruning method is set to {}, but no outgroup has been \
specified".format(args.prune))

    if args.min_support:
        if args.min_support < 0 or args.min_support > 100:
            _error("minimum support value ('--min-support') has to be \
either in percentage (1-100) or in decimal format between 0.0 - 1.0")
        elif args.min_support > 1:
            # convert from percentage to floating point
            args.min_support = args.min_support / 100

    if args.min_otu_occupancy:
        if args.min_otu_occupancy < 0 or args.min_otu_occupancy > 100:
            _error("minimum OTU occupancy ('--min-otu-occupancy') has to be \
either in percentage (1-100) or in decimal format between 0.0 - 1.0")
        elif args.min_otu_occupancy > 1:
            # convert from percentage to floating point
            args.min_otu_occupancy = args.min_otu_occupancy / 100

    if args.min_gene_occupancy:
        if args.min_gene_occupancy < 0 or args.min_gene_occupancy > 100:
            _error("minimum gene occupancy ('--min-gene-occupancy') has to be \
either in percentage (1-100) or in decimal format between 0.0 - 1.0")
        elif args.min_gene_occupancy > 1:
            # convert from percentage to floating point
            args.min_gene_occupancy = args.min_gene_occupancy / 100

    if args.trim_divergent:
        if args.trim_divergent < 0 or args.trim_divergent > 100:
            _error("the divergence threshold ('--trim-divergent') has to be \
either in percentage (1-100) or in decimal format between 0.0 - 1.0")
        elif args.trim_divergent > 1:
            # convert from percentage to floating point
            args.trim_divergent = args.trim_divergent / 100

    if args.min_len and args.min_len < 1:
        _error("minimum sequence length ('--min-len') has to be a positive \
integer (1, 2, 3, 4, ...)")

    if args.threads and args.threads < 1:
        _error("threads ('--threads') has to be a positive integer \
(1, 2, 3, 4, ...)")

    if args.min_taxa and args.min_taxa < 1:
        _error("minimum number of taxa ('--min-taxa') has to be a positive \
integer (1, 2, 3, 4, ...)")

    if args.trim_lb and args.trim_lb <= 0:
        _error("the factor for removing long branches ('--trim-lb') \
has to be a positive number")

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

    # remove short sequences
    if settings.min_len:
        log.trimmed_seqs = filtering.trim_short_seqs(msa, tree, settings.min_len)

    if filtering.too_few_otus(tree, settings.min_taxa):
        return log

    # trim long branches
    if settings.trim_lb:
        log.lbs_removed = list(filtering.prune_long_branches(tree, settings.trim_lb))

    if filtering.too_few_otus(tree, settings.min_taxa):
        return log

    # trim zero length branches
    if settings.trim_zero_len:
        filtering.trim_zero_len_branches(tree, settings.trim_zero_len)

    if filtering.too_few_otus(tree, settings.min_taxa):
        return log

    # collapse weakly supported nodes into polytomies
    if settings.min_support:
        log.collapsed_nodes = filtering.collapse_nodes(tree, settings.min_support)

    # mask monophyletic groups
    if settings.mask:
        if settings.mask == "pdist":
            tree, masked_seqs = mask_monophylies.pairwise_distance(tree)
        elif settings.mask == "longest":
            tree, masked_seqs = mask_monophylies.longest_isoform(msa, tree)
        log.monophylies_masked = masked_seqs

    # trim divergent sequences
    if settings.trim_divergent:
        log.divergent = decontamination.trim_divergent(
            tree, settings.trim_divergent, settings.include)

    if filtering.too_few_otus(tree, settings.min_taxa):
        return log

    # root by outgroup
    rooted = False # True if outgroup rooting is successful
    if settings.outgroup:
        if not settings.prune == "MO":
            tree, rooted = root.outgroup(tree, settings.outgroup)

    # root the tree by midpoint rooting
    if not rooted and settings.root:
        if settings.root == "midpoint":
            tree = root.midpoint(tree)

    # exclude taxa within the list settings.exclude
    if settings.exclude:
        tree = filtering.exclude(tree, settings.exclude)

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

    # get a list of paralogs
    log.paralogs = tree.paralogs()

    # prune paralogs
    log.orthologs = prune_paralogs(settings.prune, tree,
                                   settings.min_taxa, settings.outgroup)

    if settings.force_inclusion:
        log.orthologs = filtering.force_inclusion(
            log.orthologs, settings.force_inclusion)

    return log

def _get_orthologs(settings, directory="", dir_out=None):
    fasta_path = "{}{}".format(directory, settings.fasta_file)
    nw_path = "{}{}".format(directory, settings.nw_file)
    msa = fasta.read(fasta_path)
    nw_file = newick.read(nw_path)
    log = _run(settings, msa, nw_file)

    log.get_msas_out(dir_out)
    log.report(dir_out)
    return log

def _run_for_file_pairs(corr_files, settings, dir_in, dir_out):
    settings.fasta_file, settings.nw_file = corr_files
    return _get_orthologs(settings, dir_in, dir_out)

def file_pairs_from_directory(directory):
    """Takes the path to a directory as an input. Finds corresponding files, or
    file pairs, where each pair consists of the path to a multiple sequence
    alignment (MSA) file and a Newick tree file. These pairs are found by
    looking for files with the same name, but with different extensions, within
    the provided directory. Fasta and Newick files must have a filetype
    extension that is recognized by PhyloPyPruner (see the FASTA_EXTENSIONS and
    NW_EXTENSIONS variables for recognized extensions).

    Parameters
    ----------
    directory : string
        Path to a directory containing one or more file pairs, where each pair
        corresponds to a FASTA or Newick file.

    Returns
    -------
    corr_files : dict
        Short for 'corresponding files', that is, a file pair. This dictionary
        contains tuples of two items where the first item corresponds to the
        path to a multiple sequence alignment (MSA) file and the second item
        corresponds to the path to a Newick tree file. These pairs where
        generated by looking for pairs of files with matching filenames, but
        with different filetype extensions.
    """
    if not os.path.isdir(directory):
        _error("input directory {} does not exist".format(directory))

    # corresponding files; filename is key, values are tuple pairs where
    # MSA comes first and tree second
    corr_files = dict()
    msa_list=[]
    tree_list=[]

    for file in os.listdir(directory):
        filename, extension = os.path.splitext(file)
        # filename = ".".join(filename.split(".")[:-1])
        extension = extension.lower()

        if extension in FASTA_EXTENSIONS:
            msa_list.append(file)
        if extension in NW_EXTENSIONS:
            tree_list.append(file)
        # SORT and CHECK SIZES
    if len(msa_list) == len(tree_list):
        msa_list.sort()
        tree_list.sort()
        return list(zip(msa_list,tree_list))
    else:
        _error("mismatch in file pairs between trees ({}) and msa({}) in directory".format(len(tree_list),len(msa_list))
               )
        return None

    #     if extension in FASTA_EXTENSIONS or extension in NW_EXTENSIONS:
    #         if filename in corr_files.keys():
    #             if extension in NW_EXTENSIONS:
    #                 corr_files[filename] = corr_files[filename] + (file,)
    #             else:
    #                 corr_files[filename] = (file,) + corr_files[filename]
    #         else:
    #             corr_files[filename] = (file,)
    #
    # # remove file "pairs" where only 1 file was recovered
    # pairs_to_remove = list()
    # for filename in corr_files:
    #     if len(corr_files[filename]) < 2:
    #         pairs_to_remove.append(filename)
    #
    # for filename in pairs_to_remove:
    #     corr_files.pop(filename)
    #
    # return corr_files

def parse_args():
    "Parse the arguments provided by the user."
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-v", "-V", "--version",
                        action="version",
                        version=str(VERSION),
                        help="display the version number and exit")
    parser.add_argument("--output",
                        metavar="<directory>",
                        type=str,
                        default=None,
                        help="set output directory; DEFAULT: same directory \
                              as the input alignments")
    parser.add_argument("--wrap",
                        metavar="<max column>",
                        default=None,
                        type=int,
                        help="wrap output sequences at column <max column>; \
                              sequence data is kept at a single line by \
                              default")
    parser.add_argument("--threads",
                        metavar="<count>",
                        default=None,
                        type=int,
                        help="limit the number of logical cores used for \
                              processing files in parallel")
    parser.add_argument("--overwrite",
                        default=False,
                        action="store_true",
                        help="overwrite an older run, if it exist within \
                              the output directory")
    parser.add_argument("--no-plot",
                        default=False,
                        action="store_true",
                        help="shorten the run time by not generating any plot")

    group = parser.add_argument_group("input files (MSA and tree or directory)")
    group.add_argument("--msa",
                       metavar="<MSA>",
                       type=str,
                       help="path to a multiple sequence alignment (MSA) in \
                             FASTA format")
    group.add_argument("--tree",
                       metavar="<tree>",
                       type=str,
                       help="path to a Newick tree file")
    group.add_argument("--dir",
                       metavar="<directory>",
                       type=str,
                       default=None,
                       help="path to a directory containing multiple MSAs \
                             and trees")

    group = parser.add_argument_group("prefiltering")
    group.add_argument("--min-taxa",
                       metavar="<threshold>",
                       type=int,
                       default=4,
                       help="minimum number of OTUs allowed in output; \
                             DEFAULT: 4")
    group.add_argument("--min-len",
                       metavar="<threshold>",
                       default=None,
                       type=int,
                       help="remove sequences shorter than <threshold>")
    group.add_argument("--min-support",
                       metavar="<threshold>",
                       default=None,
                       type=float,
                       help="collapse nodes with a support value below \
                             <threshold> into polytomies; in percent \
                             or in decimal format (0.0-1.0)")
    group.add_argument("--trim-lb",
                       default=None,
                       metavar="<factor>",
                       type=float,
                       help="remove sequences with a branch length that is \
                             <factor> times larger than the standard \
                             deviation of all branches within its tree")
    group.add_argument("--min-pdist",
                       default=None,
                       metavar="<distance>",
                       type=float,
                       help="remove pairs of sequences with a tip-to-tip \
                             distance that is less than <distance>")
    group.add_argument("--include",
                       nargs="+",
                       metavar="<OTU>",
                       default=None,
                       type=str,
                       help="include the OTUs in the specified list, even if \
                             they were deemed problematic by \
                             '--trim-freq-paralogs' or '--trim-divergent'")
    group.add_argument("--exclude",
                       nargs='+',
                       metavar="<OTU>",
                       default=None,
                       type=str,
                       help="exclude these OTUs in this run")
    group.add_argument("--force-inclusion",
                       nargs="+",
                       metavar="<OTU>",
                       default=None,
                       type=str,
                       help="do not output any orthologs where these OTUs \
                             are not present")
    group.add_argument("--min-otu-occupancy",
                       metavar="<percentage>",
                       default=None,
                       type=float,
                       help="remove OTUs with less occupancy than \
                       <percentage> from the supermatrix; \
                       in percent or in decimal format (0.0-1.0)")
    group.add_argument("--min-gene-occupancy",
                       metavar="<percentage>",
                       default=None,
                       type=float,
                       help="remove genes with less occupancy than \
                       <percentage> from the supermatrix; \
                       in percent or in decimal format (0.0-1.0)")

    group = parser.add_argument_group("tree-based orthology inference")
    group.add_argument("--outgroup",
                       nargs="+",
                       metavar="<OTU>",
                       default=None,
                       type=str,
                       help="root trees using one or more outgroup OTUs if \
                             at least one outgroup OTU is present and that \
                             the ones that are present are non-repetetive \
                             and form a monophyletic group")
    group.add_argument("--root",
                       default=None,
                       type=str,
                       choices=["midpoint"],
                       help="root trees using this method in case no \
                             outgroups were provided or if outgroup rooting \
                             fails")
    group.add_argument("--mask",
                       default="pdist",
                       type=str,
                       choices=["longest", "pdist"],
                       help="specify the method for masking monophylies; \
                             DEFAULT: pairwise distance ('pdist')")
    group.add_argument("--prune",
                       default="LS",
                       type=str,
                       choices=["LS", "MI", "MO", "RT", "1to1"],
                       help="set the paralogy pruning method; default is \
                             largest subtree ('LS')")

    group = parser.add_argument_group("decontamination")
    group.add_argument("--trim-freq-paralogs",
                       default=None,
                       metavar="<factor>",
                       type=float,
                       help="remove OTUs with a paralogy frequency (PF) \
                             larger than <factor> times the standard \
                             deviation of the PF for all OTUs")
    group.add_argument("--trim-divergent",
                       default=None,
                       metavar="<percentage>",
                       type=float,
                       help="exclude OTUs on a per-alignment basis, where \
                             the ratio between the maximum pairwise distance \
                             between sequences within the OTU and the \
                             average pairwise distance between sequences \
                             outside the OTU exceeds the user-defined \
                             percentage")
    group.add_argument("--jackknife",
                       default=False,
                       action="store_true",
                       help="leave out each OTU one by one and output \
                       statistics for each OTU left out into the summary \
                       file; use the '--exclude' flag in a subsequent run \
                       to remove taxa deemed as problematic")
    group.add_argument("--subclades",
                       default=None,
                       metavar="<subclade definition file>",
                       type=str,
                       help="specify a set of subclades and analyse\
                       their overall stability")
    return parser.parse_args(args=None if sys.argv[1:] else ['--help'])

def main():
    "Parse args, run filter and infer orthologs."
    args = parse_args()
    print(ABOUT)
    _validate_arguments(args)
    summary = Summary()

    if args.subclades:
        # replace the file path with a set of taxonomic groups
        groups = taxonomic_groups.read(args.subclades)
        args.subclades = groups

    settings = Settings(args)

    if args.threads:
        threads = args.threads
    else:
        threads = cpu_count()

    dir_out = None
    if args.output:
        # create output directory if it does not currently exist
        dir_out = args.output.rstrip("/") # get rid of trailing slash in path
    elif args.dir:
        dir_out = str(args.dir).rstrip("/")
    else:
        dir_out = os.path.dirname(str(args.msa).rstrip("/"))
    dir_out = dir_out + "/phylopypruner_output"

    # if not args.overwrite and os.path.isfile(dir_out):
    if not args.overwrite and os.path.isdir(dir_out):
        question = _warning("files from a previous run exist in the output \
directory, overwrite?")
        if not _yes_or_no(question):
            exit()

    if os.path.isdir(dir_out):
        # this removes the old 'phylopypruner_output' directory
        shutil.rmtree(dir_out)

    os.makedirs(dir_out)
    os.makedirs(dir_out + ORTHOLOGS_PATH)


    with open(dir_out + LOG_PATH, "w") as log_file:
        log_file.write(ABOUT + "\n")

    with open(dir_out + ORTHO_STATS_PATH, "w") as ortho_stats_file:
        ortho_stats_file.write(ORTHOLOG_STATS_HEADER)

    with open(dir_out + HOMOLOG_STATS_PATH, "w") as homo_stats_file:
        homo_stats_file.write(HOMOLOG_STATS_HEADER)

    settings.report(args.dir, dir_out + LOG_PATH)

    if args.msa and args.tree:
        # run for a single pair of files
        summary.logs.append(_get_orthologs(settings, "", dir_out))

    elif args.dir:
        if not args.dir[-1] == "/":
            dir_in = args.dir + "/"
        else:
            dir_in = args.dir

        if not os.path.isdir(dir_in):
            _error("input directory {} does not exist".format(dir_in))

        file_pairs = file_pairs_from_directory(dir_in)

        no_of_file_pairs = len(file_pairs)
        if no_of_file_pairs < 1:
            # no file pairs found in the provided directory
            _error("no file pairs were found in the provided directory")

        if args.threads:
            threads = args.threads
        else:
            threads = cpu_count()
        pool = Pool(processes=threads)
        part_run = partial(_run_for_file_pairs, settings=settings,
                           dir_in=dir_in, dir_out=dir_out)
        progress = None

        # For debugging purposes only (gets rid of multiprocessing).
        # for pair in file_pairs:
        #     print(pair)
        #     settings.fasta_file, settings.nw_file = pair
        #     _get_orthologs(settings, dir_in, dir_out)
        # exit()

        for index, log in enumerate(pool.imap_unordered(part_run, file_pairs), 1):
            progress = "{}==>{} processing MSAs and trees{} ({}/{} file \
pairs)".format("\033[34m", "\033[0m", "\033[0m", index, no_of_file_pairs)
            print(progress, end="\r")
            sys.stdout.flush()
            summary.logs.append(log)
        pool.terminate()

        mk_sum_out_title(dir_out)
        print("")

    if not summary:
        _error("no orthologs recovered, check filetype extensions or try \
more relaxed settings")

    paralog_freq = summary.paralogy_frequency(dir_out, args.trim_freq_paralogs,
                                              args.no_plot)
    homolog_stats = summary.homolog_report(dir_out)
    otus_to_exclude = []

    if args.trim_freq_paralogs:
        freq_paralogs = decontamination.trim_freq_paralogs(
            args.trim_freq_paralogs, paralog_freq)
        if freq_paralogs:
            otus_to_exclude += freq_paralogs

    # Generate a list of OTUs to exclude and exclude them from the summary.
    if otus_to_exclude:
        if args.include:
            otus_to_exclude = [otu for otu in otus_to_exclude if not otu in
                               args.include]
        summary, ortholog_report = decontamination.prune_by_exclusion(
            summary, otus_to_exclude, dir_out, threads, homolog_stats)

    # Get OTUs and genes to exclude based on their occupancy.
    otus_to_exclude, genes_to_exclude = summary.matrix_occupancy(
        dir_out, args.min_otu_occupancy, args.min_gene_occupancy, args.no_plot)

    if otus_to_exclude:
        summary = decontamination.exclude_otus(summary, otus_to_exclude)

    if genes_to_exclude:
        summary = decontamination.exclude_genes(summary, genes_to_exclude)

    # Remove gap-only columns from the output alignments.
    summary = summary.remove_gap_only_columns()

    if settings.taxonomic_groups:
        decontamination.score_monophyly(
            summary, settings.taxonomic_groups, dir_out)

    supermatrix = Supermatrix(dir_out)
    supermatrix.partitions_from_summary(summary, dir_out)

    # Perform taxon jackknifing.
    if args.jackknife:
        decontamination.jackknife(summary, dir_out, threads)

    ortholog_report = summary.report("output", dir_out, homolog_stats)

    print(ortholog_report)

    with open(dir_out + LOG_PATH, "a") as log_file:
        log_file.write("\n" + ortholog_report)

    summary.write_msas(args.wrap)
    run_time = "Run time: {} seconds".format(round(time.time() - START_TIME, 2))
    run_time_report = "\n{}\n{}".format("-" * len(run_time), run_time)
    print(run_time_report)

    with open(dir_out + LOG_PATH, "a") as log_file:
        log_file.write("\n" + run_time_report)

if __name__ == "__main__":
    sys.path.insert(0, os.path.abspath('..'))
    START_TIME = time.time()
    main()
