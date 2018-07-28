#!/usr/bin/env python

# Author: Felix Thalen
# Date: July 24, 2018

r"""
-------------------------------------------------------------------------------
 ___ _        _     ___      ___
| _ \ |_ _  _| |___| _ \_  _| _ \_ _ _  _ _ _  ___ _ _
|  _/ ' \ || | / _ \  _/ || |  _/ '_| || | ' \/ -_) '_|
|_| |_||_\_, |_\___/_|  \_, |_| |_|  \_,_|_||_\___|_|
         |__/           |__/

A contamination aware phylogenetic tree-based software for orthology inference.

Refer to the wiki (https://github.com/fethalen/phylopypruner/wiki) for a
tutorial and explanations of the implemented algorithms.

Example:
  ./phylopypruner.py msa.fa tree.nw --min-taxa 4 --min-seq 40
-------------------------------------------------------------------------------
"""

import argparse
import fasta
import newick

VERSION = 0.1

def parse_args():
    """ Parse the arguments provided by the user.
    """
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('msa', type=str, \
            help='path to a multiple sequence alignment in FASTA format')
    parser.add_argument('tree', type=str, \
            help='path to a Newick tree file')
    parser.add_argument('-v', '--version', action='version', \
            version=str(VERSION), help='display the version number and exit')
    parser.add_argument('--output', type=str, default='.', \
            help='path to output directory; default: current directory')
    parser.add_argument('--min-taxa', type=int, default=4, \
            help='minimum number of taxa allowed in orthologs')
    parser.add_argument('--min-seq', type=int, default=40, \
            help='minimum number of positions allowed in output sequences')
    parser.add_argument('--collapse-nodes', type=int, default=80, \
            help='collapse nodes with a support value below a treshold')
    return parser.parse_args()

def main():
    """ Parse and store a Newick file, iterate over its tips and store its name
    and branch length to a dictionary. Use that dictionary to identify branches
    that deviate more than 'factor' standard deviations from the median branch
    length and output their name and distance, separated by tab.
    """
    args = parse_args()

    msa = fasta.read(args.msa)
    tree = newick.read(args.tree)

    for sequence in msa.sequences:
        print(sequence)

    # count = 0
    # for node in prune_paralogs.maximum_inclusion(tree, args.min_taxa):
    #     print('Node no. {}'.format(count))
    #     count += 1
    #     for leaf in node.iter_leaves():
    #         otu = re.split('\||@', leaf.name)[0]
    #         print(otu)

if __name__ == '__main__':
    main()
