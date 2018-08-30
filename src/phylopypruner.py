#!/usr/bin/env python
#pylint: disable=too-many-branches

# Author: Felix Thalen
# Date: July 24, 2018

r"""
------------------------------------------------------------------------------
PhyloPyPruner

Contamination aware and modular phylogenetic tree-based orthology inference.

Refer to the wiki (https://github.com/fethalen/phylopypruner/wiki) for a
tutorial and explanations of the implemented algorithms.
------------------------------------------------------------------------------
example usage:
  ./phylopypruner.py 16s.fa 16s.tre --min-taxa 4 --min-seq 200 --support 0.7
  --trim-lb 5 --outgroup Drosophila --root midpoint --mask pdist --prune MI
"""

import argparse
import os
import fasta
import newick
import filtering
import mask_monophylies
import root
from prune_paralogs import prune_paralogs
from msa import MultipleSequenceAlignment
from log import Log

VERSION = 0.1

def _validate_input(msa, tree, tree_path):
    "Test to see if MSA and tree entries matches."
    descriptions = list(msa.iter_descriptions())
    names = list(tree.iter_names())

    if set(descriptions).intersection(names) < set(descriptions):
        raise AssertionError("MSA {} does not match tree \
                {}".format(msa.filename, tree_path))

def _validate_arguments(args):
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

def _file_out(path, directory=None, index=None):
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

    if index:
        index = "_{}".format(index)

    return "{}/{}_pruned{}{}".format(directory, basename, index, extension)

def parse_args():
    """
    Parse the arguments provided by the user.
    """
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("msa",
                        type=str,
                        help="multiple sequence alignment (MSA) in FASTA \
                              format")
    parser.add_argument("tree",
                        type=str,
                        help="Newick tree file")
    parser.add_argument("-v", "--version",
                        action="version",
                        version=str(VERSION),
                        help="display the version number and exit")
    parser.add_argument("--output",
                        metavar="directory",
                        type=str,
                        default=".",
                        help="path to output directory")
    parser.add_argument("--min-taxa",
                        metavar="treshold",
                        type=int,
                        default=4,
                        help="minimum number of OTUs allowed in output")
    parser.add_argument("--min-seq",
                        metavar="treshold",
                        default=40,
                        type=int,
                        help="remove sequences shorter than 'treshold'")
    parser.add_argument("--min-support",
                        metavar="treshold",
                        default=None,
                        type=float,
                        help="collapse nodes with a support value below \
                              'treshold' into polytomies")
    parser.add_argument("--trim-lb",
                        default=None,
                        metavar="factor",
                        type=int,
                        help="remove sequences with a branch length that is \
                              'factor' times the standard deviation of all \
                              branches.")
    parser.add_argument("--outgroup",
                        nargs='+',
                        metavar="OTU",
                        default=None,
                        type=str,
                        help="one or more outgroup OTUs")
    parser.add_argument("--root",
                        default=None,
                        type=str,
                        choices=["midpoint", "clock", "MAD"],
                        help="reroot tree using this rooting method")
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
    parser.add_argument("--verbose",
                        default=False,
                        action="store_true",
                        help="show a more detailed report")
    parser.add_argument("--quiet",
                        default=False,
                        action="store_true",
                        help="don't output a report")
    return parser.parse_args()

def main():
    "Parse arguments, run filter and infer orthologs."
    args = parse_args()

    _validate_arguments(args)

    # parse the MSA and tree and create a log to keep track of operations
    nw_file = args.tree
    msa_file = args.msa
    msa = fasta.read(msa_file)
    tree = newick.read(nw_file)
    log = Log(VERSION, msa, tree, args)

    # test to see if the entries in the MSA and tree matches
    _validate_input(msa, tree, args.tree)

    # remove short sequences
    if args.min_seq:
        log.trimmed_seqs = filtering.trim_short_seqs(msa, tree, args.min_seq)

    # prune long branches
    if args.trim_lb:
        log.lbs_removed = list(filtering.prune_long_branches(tree,
                                                             args.trim_lb))

    # mask monophylies
    if args.mask:
        if args.mask == "pdist":
            tree, masked_seqs = mask_monophylies.pairwise_distance(tree)
        elif args.mask == "longest":
            tree, masked_seqs = mask_monophylies.longest_isoform(msa, tree)
        log.monophylies_masked = masked_seqs

    # collapse weakly supported nodes into polytomies
    if args.min_support:
        log.pruned_sequences = filtering.collapse_nodes(tree, args.min_support)

    # We run masked monophylies again, since some sequences might be missed
    # otherwise. Only running masked monophylies here doesn't solve the
    # problem.
    if args.mask:
        if args.mask == "pdist":
            tree, masked_seqs = mask_monophylies.pairwise_distance(tree)
        elif args.mask == "longest":
            tree, masked_seqs = mask_monophylies.longest_isoform(msa, tree)
        log.monophylies_masked.append(masked_seqs)

    # exit if number of OTUs < treshold
    if filtering.too_few_otus(tree, args.min_taxa):
        print("too few OTUs in tree {}, exiting".format(nw_file))
        exit()

    if args.outgroup:
        if not args.prune == "MO":
            rerooted_tree = root.outgroup(tree, args.outgroup)
            if rerooted_tree:
                tree = rerooted_tree

    # get a list of paralogs
    log.paralogs = tree.paralogs()

    # prune paralogs
    log.orthologs = prune_paralogs(
        args.prune, tree, args.min_taxa, args.outgroup)

    log.report(args.verbose)
    print("\n" + _file_out(msa_file))

if __name__ == "__main__":
    main()
