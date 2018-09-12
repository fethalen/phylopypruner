#!/usr/bin/env python

# Author: Felix Thalen
# Date: July 18, 2018

"""
  Generate statistics for the Newick tree files in the directory contained
  within the provided <path>.
"""

import os
import argparse
import newick

NW_EXTENSIONS = (".tre", ".nw", ".newick")

def get_nw_stats(path):
    """ Takes the path to a file in Newick format as an input and get the
    average support and
    """
    tree = newick.read(path)

    for node in tree.traverse():
        if node.is_leaf():
            leaf_count += 1
            total_dist += node.dist
        else:
            branch_count += 1
            total_support += node.support
    return leaf_count, total_dist, branch_count, total_support

def parse_args():
    """ Parse the arguments provided by the user.
    """
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("path", type=str, \
            help="The file containing the tree in Newick format")
    return parser.parse_args()

def main():
    """ Parse and store a Newick file, iterate over its tips and store its name
    and branch length to a dictionary. Use that dictionary to identify branches
    that deviate more than 'factor' standard deviations from the median branch
    length and output their name and distance, separated by tab.
    """
    args = parse_args()
    dir_path = args.path
    total_dist = 0
    total_support = 0
    leaf_count = 0
    branch_count = 0

    for file in os.listdir(dir_path):
        if file.endswith(NW_EXTENSIONS):
            nw_file = os.path.join(dir_path, file)
            get_nw_stats(nw_file)

    avg_dist = total_dist / leaf_count
    avg_support = total_support / branch_count

    print("Average distance: {}".format(avg_dist))
    print("Average support: {}".format(avg_support))

if __name__ == '__main__':
    main()
