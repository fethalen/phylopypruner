#!/usr/bin/env python

"""
Implements various paralog pruning algorithms, used for finding orthologs in a
gene tree.
"""

def _repetetive_otus(node):
    """
    Takes a TreeNode object as an input, returns true if one or more OTUs are
    present in the tree more than once.
    """
    otus = list(node.iter_otus())
    return len(otus) > len(set(otus))

def largest_subtree(tree, min_taxa):
    """
    Takes a TreeNode object and a minimum number of taxa as an input. Find and
    return the TreeNode object of the largest subtree with non-repeating OTUs.
    """
    max_subtree = None
    max_size = 0

    for node in tree.traverse_preorder():
        if not _repetetive_otus(node):
            if node.is_leaf():
                continue

            leaf_count = sum(1 for _ in node.iter_leaves())

            if leaf_count >= min_taxa:
                if not max_subtree:
                    max_subtree = node
                    max_size = leaf_count
                elif leaf_count > max_size:
                    max_subtree = node
                    max_size = leaf_count

    return max_subtree

def maximum_inclusion(tree, min_taxa):
    """
    Takes a TreeNode object and a minimum number of taxa as an input. Find and
    yield the largest subtree with non-repeating OTUs, cut it off and iterate
    over the remaining tree.
    """
    max_subtree = largest_subtree(tree, min_taxa)

    if max_subtree.is_root():
        # the entire tree consists of non-repetetive taxa only
        yield max_subtree
        max_subtree = None

    while max_subtree:
        yield max_subtree
        tree.remove_node(max_subtree)
        max_subtree = largest_subtree(tree, min_taxa)

def rooted_tree(tree, min_taxa, outgroup):
    """
    """
    for node in tree.traverse_preorder():
        if node.outgroups_only:
            pass

def monophyletic_outgroups(tree, min_taxa, outgroup):
    """
    """
    pass

def one_to_one_orthologs(tree):
    """
    Takes a Newick tree as an input and returns an unmodified tree if and only
    if the OTUs are non-repetetive.
    """
    if not _repetetive_otus(tree):
        return tree

def get_paralogs(tree):
    "Return a list of paralogs present in the provided TreeNode object."
    if not _repetetive_otus(tree):
        return

    seen = set()
    paralogs = set()

    for otu in tree.iter_otus():
        if otu in seen:
            paralogs.add(otu)
        else:
            seen.add(otu)
    return paralogs

def prune_paralogs(method, tree, min_taxa, outgroup):
    """
    Takes the name of the paralogy pruning algorithm to use and a TreeNode
    object as an input. Returns a list of TreeNode objects, where each object
    is an inferred orthology tree.
    """
    subtrees = []
    methods = ("LS", "MI", "MO", "RT", "1to1")

    if not method in methods:
        AssertionError("unknown paralogy pruning method: '{}'".format(method))

    if method == "LS":
        subtrees.append(largest_subtree(tree, min_taxa))
    elif method == "MI":
        subtrees = list(maximum_inclusion(tree, min_taxa))
    elif method == "MO":
        subtrees = list(monophyletic_outgroups(tree, min_taxa, outgroup))
    elif method == "RT":
        subtrees = list(rooted_tree(tree, min_taxa, outgroup))
    elif method == "1to1":
        if one_to_one_orthologs(tree):
            subtrees.append(tree)

    return subtrees
