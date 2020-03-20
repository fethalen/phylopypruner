#!/usr/bin/env python

# Author: Felix Thalen
# Date: July 24, 2018

"""
A tree-based orthology inference program with additional functionality for
reducing contamination. See gitlab.com/fethalen/phylopypruner for details.
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
from phylopypruner import report
from phylopypruner import taxonomic_groups
from phylopypruner.prune_paralogs import prune_paralogs
from phylopypruner.summary import Summary
from phylopypruner.summary import mk_sum_out_title
from phylopypruner.log import Log
from phylopypruner.settings import Settings
from phylopypruner.supermatrix import Supermatrix

VERSION = pkg_resources.require("phylopypruner")[0].version
FASTA_EXTENSIONS = {".fa", ".fas", ".fasta", ".fna", ".faa", ".fsa", ".ffn",
                    ".frn"}
NW_EXTENSIONS = {".newick", ".nw", ".tre", ".tree", ".out", ".treefile"}
ORTHO_STATS_PATH = "/output_alignment_stats.csv"
HOMOLOG_STATS_PATH = "/input_alignment_stats.csv"
ORTHOLOG_STATS_HEADER = "filename,otus,sequences,meanSeqLen,shortestSeq,longestSeq,\
pctMissingData,alignmentLen\n"
HOMOLOG_STATS_HEADER = "filename,otus,sequences,meanSeqLen,shortestSeq,longestSeq,\
pctMissingData,alignmentLen,shortSequencesRemoved,longBranchesRemoved,\
monophyliesMasked,nodesCollapsed,divergentOtusRemoved\n"
ORTHO_STATS_PATH = "/output_alignment_stats.csv"
LOG_PATH = "/phylopypruner.log"
ORTHOLOGS_PATH = "/output_alignments"
TIMESTAMP = datetime.datetime.now().strftime("%A, %d. %B %Y %I:%M%p")
ABOUT = report.underline("PhyloPyPruner version {}".format(VERSION))
ABOUT_LOG = "PhyloPyPruner version {}\n{}\n{}".format(
    VERSION, TIMESTAMP, "-" * len(TIMESTAMP))
NO_FILES = """You did not specify any input data. Use the flag '--dir', followed
by the path to a directory, to point to a directory which contain your
multiple sequence alignments (MSAs) and input trees."""


def _validate_input(msa, tree, tree_path):
    "Test to see if MSA and tree entries matches."
    descriptions = list(msa.iter_descriptions())
    names = list(tree.iter_names())

    if set(descriptions).intersection(names) < set(descriptions):
        print("example tree names:", names[:2], file=sys.stderr)
        print("example sequences:", descriptions[:2], file=sys.stderr)
        report.error("MSA names don't match tree \n   {}\n   {}".format(
            msa.filename, tree_path))


def _no_files(args):
    "Returns True if the required files are not provided."
    return not args.dir


def _validate_arguments(args):
    "Performs a series of checks to validate the input before execution."

    errors = False

    if _no_files(args):
        report.error(NO_FILES)
        errors = True

    if not args.outgroup and args.prune == "MO" or \
       not args.outgroup and args.prune == "RT":
        report.error("pruning method is set to {}, but no outgroup has been \
specified".format(args.prune))
        errors = True

    if args.min_support:
        if args.min_support < 0 or args.min_support > 100:
            report.error("minimum support value ('--min-support') has to be \
either in percentage (1-100) or in decimal format between 0.0 - 1.0")
            errors = True
        elif args.min_support > 1:
            # convert from percentage to floating point
            args.min_support = args.min_support / 100

    if args.min_otu_occupancy:
        if args.min_otu_occupancy < 0 or args.min_otu_occupancy > 100:
            report.error("minimum OTU occupancy ('--min-otu-occupancy') has to be \
either in percentage (1-100) or in decimal format between 0.0 - 1.0")
            errors = True
        elif args.min_otu_occupancy > 1:
            # convert from percentage to floating point
            args.min_otu_occupancy = args.min_otu_occupancy / 100

    if args.min_gene_occupancy:
        if args.min_gene_occupancy < 0 or args.min_gene_occupancy > 100:
            report.error("minimum gene occupancy ('--min-gene-occupancy') has to be \
either in percentage (1-100) or in decimal format between 0.0 - 1.0")
            errors = True
        elif args.min_gene_occupancy > 1:
            # convert from percentage to floating point
            args.min_gene_occupancy = args.min_gene_occupancy / 100

    if args.trim_divergent:
        if isinstance(args.trim_divergent, int):
            args.trim_divergent = args.trim_divergent / 100
        if args.trim_divergent < 0:
            report.error("the divergence threshold ('--trim-divergent') has to be \
a positive number")
            errors = True

    if args.min_len and args.min_len < 1:
        report.error("minimum sequence length ('--min-len') has to be a positive \
integer (1, 2, 3, 4, ...)")
        errors = True

    if args.threads and args.threads < 1:
        report.error("threads ('--threads') has to be a positive integer \
(1, 2, 3, 4, ...)")
        errors = True

    if args.min_taxa and args.min_taxa < 1:
        report.error("minimum number of taxa ('--min-taxa') has to be a positive \
integer (1, 2, 3, 4, ...)")
        errors = True

    if args.trim_lb and args.trim_lb <= 0:
        report.error("the factor for removing long branches ('--trim-lb') \
has to be a positive number")
        errors = True

    if errors:
        sys.exit()


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

    # exclude taxa within the list settings.exclude
    if settings.exclude:
        tree = filtering.exclude(tree, settings.exclude)

    if filtering.too_few_otus(tree, settings.min_taxa):
        return log

    # remove short sequences
    if settings.min_len:
        log.trimmed_seqs = filtering.trim_short_seqs(
            msa, tree, settings.min_len)

    if filtering.too_few_otus(tree, settings.min_taxa):
        return log

    # trim long branches
    if settings.trim_lb:
        log.lbs_removed = list(filtering.prune_long_branches(
            tree, settings.trim_lb))

    if filtering.too_few_otus(tree, settings.min_taxa):
        return log

    # trim zero length branches
    if settings.trim_zero_len:
        log.ultra_short_branches = filtering.trim_zero_len_branches(
            tree, settings.trim_zero_len)

    if filtering.too_few_otus(tree, settings.min_taxa):
        return log

    # collapse weakly supported nodes into polytomies
    if settings.min_support:
        log.collapsed_nodes = filtering.collapse_nodes(
            tree, settings.min_support)

    # mask monophyletic groups
    if settings.mask:
        if settings.mask == "pdist":
            tree, masked_seqs = mask_monophylies.pairwise_distance(tree)
        elif settings.mask == "longest":
            tree, masked_seqs = mask_monophylies.longest_isoform(msa, tree)
        log.monophylies_masked = masked_seqs

    # trim divergent sequences
    if settings.trim_divergent:
        log.divergent, log.divergent_removed = decontamination.trim_divergent(
            tree, settings.trim_divergent, settings.include)

    if filtering.too_few_otus(tree, settings.min_taxa):
        return log

    # root by outgroup
    rooted = False   # True if outgroup rooting is successful
    if settings.outgroup:
        if not settings.prune == "MO":
            tree, rooted = root.outgroup(tree, settings.outgroup)

    # root the tree by midpoint rooting
    if not rooted and settings.root:
        if settings.root == "midpoint":
            tree = root.midpoint(tree)

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
    file_pairs : dict
        This dictionary contains tuples of two items where the first item
        corresponds to the path to a multiple sequence alignment (MSA) file and
        the second item corresponds to the path to a Newick tree file. These
        pairs where generated by looking for pairs of files with matching
        filenames, but with different filetype extensions.
    """
    if not os.path.isdir(directory):
        report.error("input directory {} does not exist".format(directory))

    report.progress_bar("fetching files from directory")
    print("", file=sys.stderr)

    # file_pairs; filename is key, values are tuple pairs where
    # MSA comes first and tree second
    file_pairs = dict()

    for file in os.listdir(directory):
        extension = os.path.splitext(file)[1]
        extension = extension.lower()
        filename = file.split(os.extsep)[0]

        if extension in FASTA_EXTENSIONS or extension in NW_EXTENSIONS:
            if filename in file_pairs.keys():
                if extension in NW_EXTENSIONS:
                    file_pairs[filename] = file_pairs[filename] + (file,)
                else:
                    file_pairs[filename] = (file,) + file_pairs[filename]
            else:
                file_pairs[filename] = (file,)

    # Report errors if more than 2 file pairs are encountered and ignore file
    # pairs with fewer than 2 files.
    errors = False
    to_remove = list()
    for file_pair in file_pairs:
        if len(file_pairs[file_pair]) > 2:
            files = list(file_pairs[file_pair])
            message = "{} files found with the name '{}' (expected 2):".format(
                len(file_pairs[file_pair]), file_pair)
            report.error(message)
            print("  " + ", ".join([pair for pair in files]), file=sys.stderr)
            errors = True
        elif len(file_pairs[file_pair]) < 2:
            to_remove.append(file_pair)

    for file_pair in to_remove:
        file_pairs.pop(file_pair)

    if errors:
        sys.exit()

    return [file_pairs[file_pair] for file_pair in file_pairs]


def parse_args():
    "Parse the arguments provided by the user."
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-v", "-V", "--version",
                        action="version",
                        version=str(VERSION),
                        help="display the version number and exit")
    parser.add_argument("--overwrite",
                        default=False,
                        action="store_true",
                        help="overwrite pre-existing output files without\
                              asking")
    parser.add_argument("--no-plot",
                        default=False,
                        action="store_true",
                        help="do not generate any plots (faster)")
    parser.add_argument("--wrap",
                        metavar="<number>",
                        default=None,
                        type=int,
                        help="wrap output sequences at column <number>, instead\
                              of writing each sequence to a single line")
    parser.add_argument("--threads",
                        metavar="<number>",
                        default=None,
                        type=int,
                        help="use <number> threads instead of up to 10 threads")
    parser.add_argument("--output",
                        metavar="<directory>",
                        type=str,
                        default=None,
                        help="save output files to <directory>, instead of\
                              the input directory")

    group = parser.add_argument_group("input data")
    group.add_argument("--dir",
                       metavar="<directory>",
                       type=str,
                       default=None,
                       help="a <directory> containing 1+ alignment and tree files")

    group = parser.add_argument_group("prefilters")
    group.add_argument("--exclude",
                       nargs='+',
                       metavar="OTU",
                       default=None,
                       type=str,
                       help="exclude these OTUs")
    group.add_argument("--min-len",
                       metavar="<number>",
                       default=None,
                       type=int,
                       help="remove sequences which are shorter than <number> bases")
    group.add_argument("--trim-lb",
                       default=None,
                       metavar="<factor>",
                       type=float,
                       help="remove branches longer than <factor> standard\
                             deviations of all branches")
    group.add_argument("--min-pdist",
                       default=None,
                       metavar="<distance>",
                       type=float,
                       help="remove sequence pairs with less tip-to-tip\
                             distance than <distance>")
    group.add_argument("--min-support",
                       metavar="<percentage>",
                       default=None,
                       type=float,
                       help="collapse nodes with less support than\
                             <percentage> into polytomies")
    group.add_argument("--mask",
                       default="pdist",
                       type=str,
                       choices=["longest", "pdist"],
                       help="if 2+ sequences from a single OTU forms a clade,\
                             choose which sequence to keep using this method")
    group.add_argument("--trim-divergent",
                       default=None,
                       metavar="<percentage>",
                       type=float,
                       help="for each alignment: discard all sequences from an OTU on a\
                       per-alignment-basis, if the ratio between the largest\
                       pairwise distance of sequences from this OTU and\
                       the average pairwise distance of sequences from this\
                       OTU to other's exceed this <percentage>")
    group.add_argument("--trim-freq-paralogs",
                       default=None,
                       metavar="<factor>",
                       type=float,
                       help="exclude OTUs with more paralogy frequency (PF)\
                             than <factor> standard deviations of all PFs")
    group.add_argument("--include",
                       nargs="+",
                       metavar="OTU",
                       default=None,
                       type=str,
                       help="include these OTUs, even if deemed problematic\
                             by '--trim-freq-paralogs' or '--trim-divergent'")

    group = parser.add_argument_group("rooting")
    group.add_argument("--outgroup",
                       nargs="+",
                       metavar="OTU",
                       default=None,
                       type=str,
                       help="root trees using these OTUs if at least one OTU\
                             is present and if all present OTUs are\
                             non-repetetive and form a clade")
    group.add_argument("--root",
                       default=None,
                       type=str,
                       choices=["midpoint"],
                       help="root trees using this method when outgroup\
                             rooting was not performed")

    group = parser.add_argument_group("paralogy pruning")
    group.add_argument("--prune",
                       default="LS",
                       type=str,
                       choices=["LS", "MI", "MO", "RT", "1to1"],
                       help="select the paralogy pruning method (default: LS)")

    group = parser.add_argument_group("postfilters")
    group.add_argument("--force-inclusion",
                       nargs="+",
                       metavar="OTU",
                       default=None,
                       type=str,
                       help="discard output alignments where these OTUs are\
                             missing")
    group.add_argument("--min-taxa",
                       metavar="<number>",
                       type=int,
                       default=4,
                       help="discard output alignments with fewer OTUs than\
                             <number> (4 by default)")
    group.add_argument("--min-otu-occupancy",
                       metavar="<percentage>",
                       default=None,
                       type=float,
                       help="do not include OTUs with less occupancy than\
                             <percentage>")
    group.add_argument("--min-gene-occupancy",
                       metavar="<percentage>",
                       default=None,
                       type=float,
                       help="discard output alignments with less occupancy\
                             than <percentage>")

    group = parser.add_argument_group("post-processing")
    group.add_argument("--subclades",
                       default=None,
                       metavar="FILE",
                       type=str,
                       help="specify a set of subclades within this file\
                             and analyse their overall stability")
    group.add_argument("--jackknife",
                       default=False,
                       action="store_true",
                       help="exclude each OTU one by one, rerun the whole\
                             analysis and generate statistics for each\
                             subsample")
    return parser.parse_args(args=None if sys.argv[1:] else ['--help'])


def main():
    "Parse args, run filter and infer orthologs."
    args = parse_args()
    print(ABOUT, file=sys.stderr)
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
        dir_out = args.output.rstrip("/")  # get rid of trailing slash in path
    elif args.dir:
        dir_out = str(args.dir).rstrip("/")
    dir_out = dir_out + "/phylopypruner_output"

    # if not args.overwrite and os.path.isfile(dir_out):
    if not args.overwrite and os.path.isdir(dir_out):
        if not report.yes_or_no(report.warning("found files from a previous \
run, overwrite these files?", display=False)):
            exit()

    if os.path.isdir(dir_out):
        # this removes the old 'phylopypruner_output' directory
        shutil.rmtree(dir_out)

    os.makedirs(dir_out)
    os.makedirs(dir_out + ORTHOLOGS_PATH)

    with open(dir_out + LOG_PATH, "w") as log_file:
        log_file.write(ABOUT_LOG + "\n")

    with open(dir_out + ORTHO_STATS_PATH, "w") as ortho_stats_file:
        ortho_stats_file.write(ORTHOLOG_STATS_HEADER)

    with open(dir_out + HOMOLOG_STATS_PATH, "w") as homo_stats_file:
        homo_stats_file.write(HOMOLOG_STATS_HEADER)

    if not args.dir[-1] == "/":
        dir_in = args.dir + "/"
    else:
        dir_in = args.dir

    if not os.path.isdir(dir_in):
        shutil.rmtree(dir_out)
        report.error("input directory {} does not exist".format(dir_in))
        sys.exit()

    # count and report the number of threads used
    thread_count = cpu_count()
    if args.threads:
        threads = args.threads if threads <= thread_count else thread_count
    else:
        threads = thread_count if thread_count <= 10 else 10
    pool = Pool(processes=threads)
    part_run = partial(_run_for_file_pairs, settings=settings,
                       dir_in=dir_in, dir_out=dir_out)
    report.progress_bar("using {} out of {} available threads".format(
        threads, thread_count))
    print("", file=sys.stderr)

    file_pairs = file_pairs_from_directory(dir_in)

    no_of_file_pairs = len(file_pairs)
    if no_of_file_pairs < 1:
        # no file pairs found in the provided directory
        report.error("no file pairs were found in the provided directory")
        shutil.rmtree(dir_out)
        sys.exit()

    for index, log in enumerate(pool.imap_unordered(part_run, file_pairs), 1):
        message = "processing MSAs and trees ({}/{})".format(index,
                                                             no_of_file_pairs)
        report.progress_bar(message)
        summary.logs.append(log)
    pool.terminate()

    mk_sum_out_title(dir_out)
    print("")

    if not summary:
        report.error("no orthologs recovered, check filetype extensions or try \
more relaxed settings")
        sys.exit()

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
            otus_to_exclude = [otu for otu in otus_to_exclude if otu not in
                               args.include]
        summary = decontamination.prune_by_exclusion(
            summary, otus_to_exclude, dir_out, threads)

    # Get OTUs and genes to exclude based on their occupancy.
    otus_to_exclude, genes_to_exclude = summary.matrix_occupancy(
        dir_out, args.min_otu_occupancy, args.min_gene_occupancy, args.no_plot)

    if genes_to_exclude:
        summary = decontamination.exclude_genes(summary, genes_to_exclude)
        summary.discarded_genes = genes_to_exclude

    if otus_to_exclude:
        summary = decontamination.exclude_otus(summary, otus_to_exclude)
        summary.discarded_otus = otus_to_exclude

    # Remove gap-only columns from the output alignments.
    summary = summary.remove_gap_only_columns()

    if settings.taxonomic_groups:
        decontamination.score_monophyly(
            summary, settings.taxonomic_groups, dir_out)

    # Perform taxon jackknifing.
    if args.jackknife:
        decontamination.jackknife(summary, dir_out, threads)

    ortholog_report = summary.report("output", dir_out, homolog_stats)
    summary.write_msas(args.wrap)

    # concatenate output alignments into a supermatrix
    supermatrix = Supermatrix(dir_out)
    supermatrix.partitions_from_summary(summary, dir_out)

    # print the output
    path_out = report.print_path(dir_out, display=False)
    report.progress_bar("wrote output to:\n  {}\n".format(path_out))

    # print settings
    print("", file=sys.stderr)
    settings.print_settings()

    # print alignment statistics
    print(ortholog_report)
    run_time = "\ncompleted in {} seconds".format(
        round(time.time() - START_TIME, 2))

    with open(dir_out + LOG_PATH, "a") as log_file:
        ortholog_report = ortholog_report.replace("\33[0m", "")
        ortholog_report = ortholog_report.replace("\33[4m", "")
        log_file.write("\n" + ortholog_report)
        log_file.write(
            "\n\nReuse these parameters:\n" +
            " ".join([arg for arg in sys.argv]))
        log_file.write("\n" + run_time)


if __name__ == "__main__":
    sys.path.insert(0, os.path.abspath('..'))
    START_TIME = time.time()
    main()
