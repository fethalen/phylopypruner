#!/usr/bin/env python
#pylint: disable=too-many-branches

"""
Implements various algorithms for monophyly masking, used to prune away
in-paralogs, isoforms, or other moonophyletic groups within the same OTU.
"""

from collections import defaultdict
from tree_node import TreeNode

def _get_new_root(node):
    # workaround to get rid of empty nodes at root
    new_root = node
    for branch in node.iter_branches():
        if branch.is_root:
            if not len(branch) is 1:
                new_root = branch
                new_root = TreeNode(branch.name, branch.dist)
                new_root.children = branch.children
                break

    leaves_to_remove = set()
    for leaf in new_root.iter_leaves():
        if not leaf.name:
            leaves_to_remove.add(leaf)
    new_root = new_root.remove_nodes(leaves_to_remove)

    while new_root.empty_leaves():
        for leaf in new_root.iter_leaves():
            if not leaf.name:
                leaf.delete()
                break

    return new_root

def _monophyletic_sister(node):
    """
    Takes a TreeNode object as an input. If both this node and the sister group
    to this node contains a non-repetetive OTU, then replace the sister group
    to this node with this node.
    """
    parent = node.parent
    if not parent:
        return

    ingroup_otu = ""

    for child in node.children:
        if child.is_leaf() and ingroup_otu:
            # non-unique OTUs found in node
            return
        elif child.is_leaf():
            ingroup_otu = child.otu()

    if not ingroup_otu:
        # no leaf within the node's children
        return

    outgroup_otu = ""
    removed = None

    for child in parent.children:
        if child.is_leaf() and outgroup_otu:
            # non-unique OTUs found in a sister group to this node
            return
        elif child.is_leaf():
            outgroup_otu = child.otu()
            if outgroup_otu == ingroup_otu:
                removed = child
                child.delete()

    return removed

def _monophyletic_polytomy(node):
    """
    Takes a TreeNode object and a string as an input. Returns tips in that node
    that contain repetetive OTUs as a dictionary, where the key is the OTU and
    the value is a list of TreeNode objects that belong to that repeated OTU.
    """
    seen = set()
    flag = set()
    multiples = defaultdict(list)

    for child in node.children:
        if not child.is_leaf():
            continue

        otu = child.otu()
        if otu in seen:
            flag.add(otu)
        seen.add(otu)

    for child in node.children:
        if not child.is_leaf():
            continue

        otu = child.otu()
        if otu in flag:
            multiples[otu].append(child)
    return multiples

def _mask(node, keep):
    """
    Takes a TreeNode object and a string as an input. If the name of a leaf
    within the branch does not match the string, then delete that leaf. Returns
    a list of removed nodes.
    """
    leaves_to_remove = set()

    for leaf in node.iter_leaves():
        if leaf.name is not keep:
            leaves_to_remove.add(leaf)
    node = node.remove_nodes(leaves_to_remove)
    return leaves_to_remove

def _sequence_len(msa, description):
    """
    Takes an MSA and a description as an input and returns the sequence length
    of that description.
    """
    sequence = msa.get_sequence(description)

    if sequence:
        return len(sequence.ungapped())

def longest_isoform(msa, node):
    """
    Takes a TreeNode object as an input. For each node that comprise a
    monophyletic group of sequences from a single OTU: compare the length of
    each sequence and remove all but the longest sequence form the tree.
    """
    removed = []
    longest_seq = 0
    keep = None

    for branch in node.iter_branches():

        if not branch.is_monophyletic():
            # branch is polyphyletic, look for polytomies with repetetive OTUs
            multiples = _monophyletic_polytomy(branch)
            for otu in multiples:
                for leaf in multiples[otu]:
                    seq_len = _sequence_len(msa, leaf.name)
                    if seq_len > longest_seq:
                        longest_seq = seq_len
                        keep = leaf.name

                leaves_to_remove = set()
                for leaf in multiples[otu]:
                    if leaf.name is not keep:
                        leaves_to_remove.add(leaf)

                node = node.remove_nodes(leaves_to_remove)
                removed = removed + list(leaves_to_remove)

        else:
            # branch is monophyletic, find the longest sequence
            for leaf in branch.iter_leaves():
                seq_len = _sequence_len(msa, leaf.name)
                if seq_len > longest_seq:
                    longest_seq = seq_len
                    keep = leaf.name

            removed += _mask(branch, keep)

        longest_seq = 0
        keep = None

    for branch in node.iter_branches():
        if branch.parent:
            monophyletic_sister = _monophyletic_sister(branch)
            if monophyletic_sister:
                removed.append(monophyletic_sister)

    new_root = _get_new_root(node)

    return new_root, removed

def pairwise_distance(node):
    """
    Takes a TreeNode object as an input. For each branch that comprise a
    monophyletic group of sequences from a single OTU: compare the distance of
    each pair of leaves within that branch and return the name of the leaf with
    the shortest distance to its sister group.
    """
    removed = []
    shortest = None

    for branch in node.iter_branches():
        if not branch.is_monophyletic():
            multiples = _monophyletic_polytomy(branch)
            for otu in multiples:
                for leaf in multiples[otu]:
                    distance = leaf.distance_to(branch)
                    if not shortest or distance < shortest:
                        keep = leaf.name
                        shortest = distance

                leaves_to_remove = set()
                for leaf in multiples[otu]:
                    if leaf.name is not keep:
                        leaves_to_remove.add(leaf)

                node = node.remove_nodes(leaves_to_remove)
                removed = removed + list(leaves_to_remove)
        else:
            # find the leaf with the shortest distance to the current branch
            for leaf in branch.iter_leaves():
                distance = leaf.distance_to(branch)

                if not shortest or distance < shortest:
                    keep = leaf.name
                    shortest = distance

            removed += _mask(branch, keep)
        keep = None

    for branch in node.iter_branches():
        if branch.parent:
            monophyletic_sister = _monophyletic_sister(branch)
            if monophyletic_sister:
                removed.append(monophyletic_sister)

    new_root = _get_new_root(node)

    return new_root, removed
