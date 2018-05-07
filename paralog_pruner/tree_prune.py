#!/usr/bin/env python

# Author: Felix Thalén
# Date: March 8, 2018

"""
  Parse the Newick tree in <path>, identify tips that are longer than
  <factor> times the standard deviations of all tips and report their name and
  branch length, separated by tab, to standard output (stdout).
"""

import argparse

def mean(data):
    """ Return the sample arithmetic mean of data, a sequence of real-valued
    numbers. The average of the empty list, '[]', is 0.
    """
    return float(sum(data)) / max(len(data), 1)

def sdm(data):
    """ Return the sum of square deviations of data.
    """
    return sum((x - mean(data))**2 for x in data)

def std(data):
    """ Return the population standard deviation of data.
    """
    if len(data) < 2:
        raise ValueError('variance requires at least two data points')
    return (sdm(data) / len(data)) ** 0.5

def find_long_branches(descendants, factor):
    """ Takes a dictionary of descendents and a treshold as an input, returns
    a list of labels, containing only the descendants that has a branch length
    that is 'factor' times longer than the standard deviation.
    """
    treshold = std(list(descendants.values())) * factor
    for descendant in descendants:
        branch_len = descendants[descendant]
        if branch_len > treshold:
            yield descendant, branch_len

def parse_args():
    """ Parse the arguments provided by the user.
    """
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('path', type=str, \
            help='The file containing the tree in Newick format')
    parser.add_argument('factor', type=float, \
            help='The number of standard deviations to use as a treshold.')
    return parser.parse_args()

def main():
    """ Parse and store a Newick file, iterate over its tips and store its name
    and branch length to a dictionary. Use that dictionary to identify branches
    that deviate more than 'factor' standard deviations from the median branch
    length and output their name and distance, separated by tab.
    """
    args = parse_args()

    import newick
    tree = newick.read(args.path)

    leaves = dict()
    for leaf in tree.iter_leaves():
        leaves[leaf.name] = leaf.dist

    for branch in find_long_branches(leaves, args.factor):
        name, dist = branch
        print("{}\t{}".format(name, dist))

if __name__ == '__main__':
    main()
