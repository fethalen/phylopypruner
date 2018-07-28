#!/usr/bin/env python

"""
Implements various algorithms for monophyly masking, used to prune away
in-paralogs, isoforms, or other moonophyletic groups within the same OTU.
"""

def longest_isoform(tree):
    """
    Takes a tree node as an input. For each node that comprise a
    monophyletic group of sequences from a single OTU: compare the length
    of each sequence and return the longest one.
    """
    for node in tree.traverse_preorder():
        if _monophyletic_otu(node):
            for leaf in node.iter_leaves:
                print(leaf)

def pairwise_distance(tree):
    """
    Takes a tree node as an input. For each node that comprise a
    monophyletic group of sequences from a single OTU: compare the length
    of each sequence and return the sequence with the shortest pairwise
    distance to its sister taxa.
    """
    pass

def _monophyletic_otu(tree):
    """
    Returns true if the provided tree only contains sequences from the same
    otu.
    """
    return len(set(tree.otus)) == 1
