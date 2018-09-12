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
  ./phylopypruner.py --msa 16s.fa --tree 16s.tre --min-taxa 4 --min-seq 200
    --min-support 0.7 --trim-lb 5 --outgroup Drosophila --root midpoint
    --mask pdist --prune MI

"""

import argparse
import os
import sys
import datetime
import fasta
import newick
import filtering
import mask_monophylies
import root
from prune_paralogs import prune_paralogs
from msa import MultipleSequenceAlignment
from log import Log
from settings import Settings

VERSION = 0.1
NO_FILES = """Please provide either a multiple sequence alignment (MSA) and a \
Newick tree,\nor a path to a directory containing multiple MSAs and Newick \
trees. Run\nPhyloPyPruner using the '-h' or '--help' flag for additional \
instructions."""
FASTA_EXTENSIONS = {".fa", ".fas", ".fasta", ".fna", ".faa", ".fsa", ".ffn",
                    ".frn"}
NW_EXTENSIONS = {".newick", ".nw", ".tre"}

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
    if _no_files(args):
        print(NO_FILES)
        exit()
    if not args.outgroup and args.prune == "MO" or \
        not args.outgroup and args.prune == "RT":
        print("No outgroup has been specified")
        exit()

def _get_sequences(msa, tree):
    """
    Takes a MultipleSequenceAlignment object and a TreeNode object as an input.
    Returns a MultipleSequenceAlignment object that contains the subset of
    sequences in the MSA that includes all sequences within the tree, provided
    that they do exist.
    """
    msa_out = MultipleSequenceAlignment()

    for name in tree.iter_names():
        match = msa.get_sequence(name)

        if not match:
            continue

        msa_out.add_sequence(match)

    return msa_out

def _file_out(path, directory=None, index=""):
    """
    Takes the path to an MSA and an optional path to a directory as an input.
    Extracts the base name and extension from the provided path. If no directory
    has been specified, then the directory is also extracted from the path.
    Returns the path to a file in the following format:
      <directory>/<basename>_pruned<extension>

    If an index has been provided, then the output will be in the following
    format:
      <directory>/<basename>_pruned_<index><extension>
    """
    if not directory:
        directory = os.path.dirname(path)
    filename = os.path.basename(path)
    basename, extension = os.path.splitext(filename)

    dir_out = directory + "/orthologs"
    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)

    if index:
        index = "_{}".format(index)

    return "{}/{}_pruned{}{}".format(dir_out, basename, index, extension)

def _run(settings, msa, tree):
    """
    Takes a dictionary that specifies a set of settings as an input. The
    dictionary should contain the following keys: msa, tree, min_taxa, min_seq,
    min_support, trim_lb, outgroup, root, mask, and prune. Runs PhyloPyPruner
    once using these settings.
    """
    # test to see if the entries in the MSA and tree matches
    _validate_input(msa, tree, settings.nw_file)

    log = Log(VERSION, msa, tree, settings)
    min_taxa = settings.min_taxa
    min_seq = settings.min_seq
    min_support = settings.min_support
    trim_lb = settings.trim_lb
    outgroup = settings.outgroup
    # rooting_method = settings.root
    pruning_method = settings.prune

    # remove short sequences
    if min_seq:
        log.trimmed_seqs = filtering.trim_short_seqs(msa, tree, min_seq)

    # trim long branches
    if trim_lb:
        log.lbs_removed = list(filtering.prune_long_branches(tree, trim_lb))

    # collapse weakly supported nodes into polytomies
    if min_support:
        log.pruned_sequences = filtering.collapse_nodes(tree, min_support)

    # mask monophyletic groups
    if settings.mask:
        if settings.mask == "pdist":
            tree, masked_seqs = mask_monophylies.pairwise_distance(tree)
        elif settings.mask == "longest":
            tree, masked_seqs = mask_monophylies.longest_isoform(msa, tree)
        log.monophylies_masked = masked_seqs

    if outgroup:
        if not pruning_method == "MO":
            tree = root.outgroup(tree, outgroup)

    # exit if number of OTUs < threshold
    if min_taxa:
        if filtering.too_few_otus(tree, min_taxa):
            print("too few OTUs in tree {}".format(settings.nw_file))
            return log

    tree.view()

    # get a list of paralogs
    log.paralogs = tree.paralogs()

    # prune paralogs
    log.orthologs = prune_paralogs(pruning_method, tree, min_taxa, outgroup)

    return log

def _get_orthologs(settings, directory="", dir_out=None, wrap=None,
                   verbose=False):
    extension_out = os.path.splitext(settings.fasta_file)[1]
    fasta_path = "{}{}".format(directory, settings.fasta_file)
    nw_path = "{}{}".format(directory, settings.nw_file)
    msa = fasta.read(fasta_path)
    nw_file = newick.read(nw_path)
    log = _run(settings, msa, nw_file)
    log.report(verbose)

    print("")
    for index, ortholog in enumerate(log.orthologs):
        if len(log.orthologs) is 1:
            file_out = _file_out(fasta_path, dir_out)
        else:
            file_out = _file_out(fasta_path, dir_out, index + 1)
        if ortholog:
            msa_out = MultipleSequenceAlignment(file_out, extension_out)
            for leaf in ortholog.iter_leaves():
                seq = msa.get_sequence(leaf.name)
                if seq:
                    msa_out.add_sequence(seq)
            if len(msa_out) > 0:
                fasta.write(msa_out, wrap)
                print("wrote: {}".format(msa_out))

def _auto():
    # configurations = itertools.product()
    pass

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
                        default=4,
                        help="minimum number of OTUs allowed in output")
    parser.add_argument("--min-seq",
                        metavar="<threshold>",
                        default=40,
                        type=int,
                        help="remove sequences shorter than <threshold>")
    parser.add_argument("--min-support",
                        metavar="<threshold>",
                        default=None,
                        type=float,
                        help="collapse nodes with a support value below\
                              <threshold> into polytomies")
    parser.add_argument("--trim-lb",
                        default=None,
                        metavar="<factor>",
                        type=int,
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
                        choices=["midpoint", "clock", "MAD"],
                        help="reroot tree using this rooting method if an\
                              outgroup hasn't been provided or if outgroup\
                              rooting fails")
    parser.add_argument("--mask",
                        default=None,
                        type=str,
                        choices=["longest", "pdist"],
                        help="mask monophylies using this method")
    parser.add_argument("--prune",
                        default=None,
                        type=str,
                        choices=["LS", "MI", "MO", "RT", "1to1"],
                        help="prune paralogs using this method")
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

    dir_out = None
    if args.output:
        # create output directory if it does not currently exist
        dir_out = args.output.rstrip("/") # get rid of trailing slash in path
        if not os.path.isdir(dir_out):
            os.makedirs(dir_out)

    if args.msa and args.tree:
        # run for a single pair of files
        print("PhyloPyPruner version {}".format(VERSION))
        print(datetime.datetime.now().strftime("%A, %d. %B %Y %I:%M%p") + "\n")
        _get_orthologs(settings, "", dir_out, args.wrap, args.verbose)
    elif args.dir:
        # run for multiple files in directory
        if not args.dir[-1] == "/":
            dir_in = args.dir + "/"
        else:
            dir_in = args.dir
        if not os.path.isdir(dir_in):
            print("input directory {} does not exist".format(dir_in))
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

        print("PhyloPyPruner version {}".format(VERSION))
        print(datetime.datetime.now().strftime("%A, %d. %B %Y %I:%M%p") + "\n")
        for pair in corr_files:
            settings.fasta_file, settings.nw_file = corr_files[pair]
            _get_orthologs(settings, dir_in, dir_out, args.wrap, args.verbose)
            print("")

if __name__ == "__main__":
    main()
