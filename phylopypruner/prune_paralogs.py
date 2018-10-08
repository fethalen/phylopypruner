#!/usr/bin/env python

"""
Implements various paralog pruning algorithms, used for finding orthologs in a
gene tree.
"""

from __future__ import absolute_import
from phylopypruner import root

def _has_enough_taxa(node, min_taxa):
    """
    Takes a TreeNode object and a treshold as an input. Returns False if there
    are fewer than [min_taxa] OTUs within the node, else True.
    """
    otus = set([_ for _ in node.iter_otus()])
    return len(otus) >= min_taxa

def largest_subtree(node, min_taxa):
    """
    Takes a TreeNode object and a minimum number of taxa as an input. Find and
    return the TreeNode object of the largest subtree with non-repeating OTUs.
    """
    max_subtree = None
    largest = 0

    for branch in node.iter_branches():
        if bool(branch.paralogs()):
            continue

        size = sum(1 for _ in branch.iter_leaves())
        if min_taxa:
            if size < min_taxa:
                continue

        if not max_subtree:
            max_subtree = branch
            largest = size
        elif size > largest:
            max_subtree = branch
            largest = size

    return max_subtree

def maximum_inclusion(tree, min_taxa):
    """
    Takes a TreeNode object and a minimum number of taxa as an input. Find and
    yield the largest subtree with non-repeating OTUs, cut it off and iterate
    over the remaining tree.
    """
    max_subtree = largest_subtree(tree, min_taxa)

    while max_subtree:
        yield max_subtree

        if max_subtree.is_root():
            break

        max_subtree.delete()

        # get rid of empty leaves that may occur after pruning
        while tree.empty_leaves():
            for leaf in tree.iter_leaves():
                if not leaf.name:
                    leaf.delete()
                    break

        max_subtree = largest_subtree(tree, min_taxa)

def rooted_tree(tree, min_taxa, outgroups):
    """
    Takes a TreeNode object, the minimum number of OTUs allowed in the output
    and a list of outgroup(s) as an input. Finds the largest subtree with the
    highest number of ingroup taxa and returns the largest subtree with
    non-repetitive OTUs.
    """
    if not tree.outgroups_present(outgroups):
        # outgroups absent
        return

    keep = None

    # find the branch with the highest number of ingroups
    for subtree in tree.iter_branches():
        if not _has_enough_taxa(subtree, min_taxa):
            continue

        if subtree.outgroups_present(outgroups):
            continue

        size = len([_ for _ in subtree.iter_leaves()])

        if not keep:
            keep = subtree
        elif size > len([_ for _ in keep.iter_leaves()]):
            keep = subtree

    max_subtree = largest_subtree(keep, min_taxa)
    return max_subtree

def monophyletic_outgroups(tree, min_taxa, outgroups):
    """
    Takes a TreeNode object, the minimum number of taxa allowed and a list of
    outgroups as an input. Looks for subtrees where all outgroup OTUs are
    present and forms a monophyletic group, cuts the tree off, roots it and
    returns the largest subtree with non-repetetive taxa within that tree. An
    iterator object is returned, where each object is the most inclusive
    subtree found within each monophyletic group.
    """
    if not tree.outgroups_present(outgroups):
        return

    for node in tree.traverse_preorder():
        # find branches where all outgroups are represented and non-repetetive
        if not node.outgroups_present(outgroups):
            continue

        if node.repetetive_outgroups(outgroups):
            continue

        if not _has_enough_taxa(node, min_taxa):
            continue

        # find the branch where all outgroups are present
        for branch in node.iter_branches():
            if not branch.outgroups_present(outgroups):
                continue

            if branch.is_monophyletic_outgroup(outgroups):
                monophyletic_outgroup = root.outgroup(branch, outgroups)
                max_subtree = largest_subtree(monophyletic_outgroup, min_taxa)
                if max_subtree:
                    yield max_subtree

def one_to_one_orthologs(tree):
    """
    Takes a Newick tree as an input and returns an unmodified tree if and only
    if the OTUs are non-repetetive.
    """
    if not bool(tree.paralogs()):
        return tree

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
        tree = largest_subtree(tree, min_taxa)
        if tree:
            subtrees.append(tree)
    elif method == "MI":
        trees = list(maximum_inclusion(tree, min_taxa))
        if trees:
            subtrees = trees
    elif method == "MO":
        trees = list(monophyletic_outgroups(tree, min_taxa, outgroup))
        if trees:
            subtrees = trees
    elif method == "RT":
        trees = rooted_tree(tree, min_taxa, outgroup)
        if trees:
            subtrees.append(trees)
    elif method == "1to1":
        if one_to_one_orthologs(tree):
            subtrees.append(tree)

    return subtrees
