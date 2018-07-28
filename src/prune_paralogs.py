#!/usr/bin/env python

"""
Implements various paralog pruning algorithms, used for finding orthologs in a
gene tree.
"""

def maximum_inclusion(tree, min_taxa):
    """
    Takes a Newick homolog tree and a minimum number of taxa treshold as an
    input. Find for the maximum inclusive subtree, cut it off and return its
    root as a putative ortholog.
    """
    for node in tree.traverse_preorder():
        if not _repetetive_otus(node):
            leaf_count = sum(1 for _ in node.iter_leaves())
            if leaf_count >= min_taxa:
                yield node

def one_to_one_orthologs(tree):
    """
    Takes a Newick tree as an input and returns an unmodified tree if and only
    if the OTUs are non-repetetive.
    """
    if not _repetetive_otus(tree):
        return tree

def _repetetive_otus(tree):
    """
    Takes a tree node as an input, returns true if one or more OTUs are
    present in the tree more than once.
    """
    return len(tree.otus) > len(set(tree.otus))
