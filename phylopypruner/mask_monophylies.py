#!/usr/bin/env python

"""
Implements various algorithms for monophyly masking, used to prune away
in-paralogs, isoforms, or other moonophyletic groups within the same OTU.
"""

from __future__ import absolute_import
from collections import defaultdict
from phylopypruner.filtering import rm_empty_root

def _monophyletic_sister(node):
    """
    Takes a TreeNode object as an input. If both this node and the sister group
    to this node contains a non-repetetive OTU, then replace the sister group
    to this node with this node.
    """
    parent = node.parent
    if not parent:
        return

    ingroup_otus = set()

    for child in node.children:
        if child.is_leaf():
            ingroup_otus.add(child.otu())

    if not ingroup_otus:
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
            if outgroup_otu in ingroup_otus:
                removed = child

    return removed

def _sequence_len(msa, description):
    """
    Takes an MSA and a description as an input and returns the sequence length
    of that description.
    """
    sequence = msa.get_sequence(description)

    if sequence:
        return len(sequence.ungapped())

def _longest_seq(node, msa):
    """
    Takes a TreeNode object and a MultipleSequenceAlignment object as an input
    and returns the node with the longest sequence within that node.
    """
    longest = None
    for leaf in node.iter_leaves():
        seq_len = _sequence_len(msa, leaf.name)
        if not longest or seq_len > longest:
            longest = seq_len
            leaf_to_keep = leaf
    return leaf_to_keep

def _smallest_distance(node):
    """
    Takes a TreeNode object as an input and returns the leaf with with the
    smallest distance to the given node.
    """
    smallest_dist = None
    for leaf in node.iter_leaves():
        dist = leaf.distance_to(node)
        if not smallest_dist or dist < smallest_dist:
            leaf_to_keep = leaf
            smallest_dist = dist
    return leaf_to_keep

def _monophyletic_branch(node, msa=None):
    no_of_leaves = len(set(node.iter_leaves()))
    keep = None

    if no_of_leaves is 2:
        # repetetive OTUs within a bifurcating node with two leaves
        if msa:
            keep = _longest_seq(node, msa)
        else:
            keep = _smallest_distance(node)
    elif no_of_leaves > 2:
        leaves_only = True
        for child in node.children:
            if not child.is_leaf():
                leaves_only = False
        if leaves_only:
            # branch is a polytomy, where each child is a leaf and
            # belongs to the same OTU
            if msa:
                keep = _longest_seq(node, msa)
            else:
                keep = _smallest_distance(node)
    if keep:
        return node.leaves_except(keep)
    else:
        return set()

def _monophyletic_polytomy(node, msa=None):
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

    shortest = None
    longest = None
    leaves_to_remove = set()

    for otu in multiples:
        for leaf in multiples[otu]:
            if msa:
                seq_len = _sequence_len(msa, leaf.name)
                if not longest or seq_len > longest:
                    keep = leaf.name
                    longest = seq_len
            else:
                distance = leaf.distance_to(node)
                if not shortest or distance < shortest:
                    keep = leaf.name
                    shortest = distance

        for leaf in multiples[otu]:
            if leaf.name is not keep:
                leaves_to_remove.add(leaf)
    return leaves_to_remove

def longest_isoform(msa, node):
    """
    Takes a TreeNode object as an input. Finds branches with repetetive OTUs
    and masks them by replacing the repetetive sequences with the longest
    sequence out of those sequences.
    """
    masked = set()
    leaves_to_remove = set()
    monophylies = True

    while monophylies:
        masked.update(leaves_to_remove)
        for leaf in leaves_to_remove:
            leaf.delete()
        # node.remove_nodes(leaves_to_remove)
        leaves_to_remove = set()
        for branch in node.iter_branches():
            child = _monophyletic_sister(branch)
            if child:
                leaves_to_remove.add(child)
            if branch.is_monophyletic():
                leaves_to_remove = _monophyletic_branch(branch, msa)
            else:
                # look for one or more repetetive OTUs within the polytomy
                leaves_to_remove.update(_monophyletic_polytomy(branch, msa))
            if leaves_to_remove:
                break
        monophylies = bool(leaves_to_remove)

    node = rm_empty_root(node)
    return node, masked

def pairwise_distance(node):
    """
    Takes a TreeNode object as an input. Finds branches with repetetive OTUs
    and masks them by replacing the repetetive sequences with the sequence with
    the shortest pairwise distance.
    """
    masked = set()
    leaves_to_remove = set()
    monophylies = True

    while monophylies:
        masked.update(leaves_to_remove)
        for leaf in leaves_to_remove:
            leaf.delete()
        node.remove_nodes(leaves_to_remove)
        leaves_to_remove = set()
        for branch in node.iter_branches():
            child = _monophyletic_sister(branch)
            if child:
                leaves_to_remove.add(child)
            if branch.is_monophyletic():
                leaves_to_remove = _monophyletic_branch(branch)
            else:
                # look for one or more repetetive OTUs within the polytomy
                leaves_to_remove.update(_monophyletic_polytomy(branch))
            if leaves_to_remove:
                break
        monophylies = bool(leaves_to_remove)

    node = rm_empty_root(node)
    return node, masked
