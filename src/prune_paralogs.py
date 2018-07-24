#!/usr/bin/env python

# Author: Felix Thalen
# Date: July 18, 2018

"""
  Parse the Newick tree in <path>.
"""

import argparse
import newick

def parse_args():
    """ Parse the arguments provided by the user.
    """
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('path', type=str, \
            help='The file containing the tree in Newick format')
    return parser.parse_args()

def main():
    """ Parse and store a Newick file, iterate over its tips and store its name
    and branch length to a dictionary. Use that dictionary to identify branches
    that deviate more than 'factor' standard deviations from the median branch
    length and output their name and distance, separated by tab.
    """
    args = parse_args()

    tree = newick.read(args.path)

    for leaf in tree.iter_leaves():
        print("{}\t{}".format(leaf.name, leaf.dist))

if __name__ == '__main__':
    main()
