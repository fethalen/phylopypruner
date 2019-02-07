"""
Tools for filtering tree nodes and sequences.
"""

from __future__ import absolute_import
import copy
from phylopypruner.tree_node import TreeNode

def _is_short_sequence(sequence, threshold):
    """
    Return true if the provided sequence is shorter than the provided threshold.

    Parameters
    ----------
    sequence : Sequence object
        The sequence you wish to consider.
    threshold : int
        Minimum number of positions allowed in sequence.

    Returns
    -------
    True of False
        True if the sequence length is shorter than the provided threshold.
    """
    return len(sequence.ungapped()) < threshold

def too_few_otus(tree, threshold):
    """
    Return true if the provided tree node object have fewer OTUs than the
    provided threshold.

    Parameters
    ----------
    tree : TreeNode object
        The tree you wish to consider.
    threshold : int
        Minimum number of OTUs allowed in the provided tree.

    Returns
    -------
    True or False
        True if the tree contains to few OTUs.
    """
    return len(set(tree.iter_otus())) < threshold

def _leaves_to_exclude(node, otus):
    if not node:
        return False
    for leaf in node.iter_leaves():
        if leaf.otu() in otus and not leaf.is_root():
            return True

def rm_empty_root(node):
    "Remove an empty node at the root of the provided node."
    if not node:
        return node

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

    return new_root.remove_nodes(leaves_to_remove)

def exclude(node, otus):
    """
    Takes a TreeNode object and a list of OTUs to exclude as an input. Returns
    the same TreeNode object but with the OTUs in the list removed.
    """
    node_excluded = copy.copy(node)
    while _leaves_to_exclude(node_excluded, otus):
        for leaf in node_excluded.iter_leaves():
            if leaf.otu() in otus and not leaf.is_root():
                leaf.delete()
                break

    if node_excluded:
        while node_excluded.empty_leaves():
            for leaf in node_excluded.iter_leaves():
                if not leaf.name:
                    leaf.delete()
                    break

    return rm_empty_root(node_excluded)

def force_inclusion(trees, otus):
    """Takes a list of TreeNode objects and a list of OTUs as an input. Returns
    the subset of these TreeNode objects where all OTUs within the list are
    present.
    """
    otus = set(otus)
    trees_w_otus = set()

    for tree in trees:
        otus_present = set()

        for leaf in tree.iter_leaves():
            otus_present.add(leaf.otu())

        if otus.issubset(otus_present):
            trees_w_otus.add(tree)

    return trees_w_otus

def _short_seqs(msa, tree, threshold):
    """
    Takes an MSA object, a TreeNode object and a threshold as an input. Returns
    True if there are sequences that are shorter than the provided threshold
    within the MSA.
    """
    for leaf in tree.iter_leaves():
        match = msa.get_sequence(leaf.name)

        if match:
            sequence = match
            if _is_short_sequence(sequence, threshold):
                return True
    return False

def trim_short_seqs(msa, tree, threshold):
    """
    Takes a TreeNode object, an MSA object and a threshold as an input. Remove
    sequences that are shorter than the provided threshold from both the MSA and
    the tree.
    """
    nodes_to_remove = set()

    while _short_seqs(msa, tree, threshold):
        for leaf in tree.iter_leaves():
            match = msa.get_sequence(leaf.name)

            if match:
                sequence = match

                if _is_short_sequence(sequence, threshold):
                    nodes_to_remove.add(leaf)

        tree.remove_nodes(nodes_to_remove)

    return nodes_to_remove

def _mean(data):
    """ Return the sample arithmetic mean of data, a sequence of real-valued
    numbers. The average of the empty list, '[]', is 0.
    """
    return float(sum(data)) / max(len(data), 1)

def _sdm(data):
    """ Return the sum of square deviations of data.
    """
    return sum((x - _mean(data))**2 for x in data)

def _std(data):
    "Return the population standard deviation of data."
    if len(data) < 2:
        raise ValueError('variance requires at least two data points')
    return (_sdm(data) / len(data)) ** 0.5

def prune_long_branches(node, factor):
    """
    Takes a TreeNode object and a integer as an input. Remove leaves that has a
    branch length that deviate more than 'factor' times the standard deviation
    of the branch length for all leaves within the nodes. Returns removed
    leaves as an iterator object.
    """
    dists = []
    leaves_to_remove = set()

    for leaf in node.iter_leaves():
        dists.append(node.distance_to(leaf))
    if len(dists) < 2:
        return leaves_to_remove
    threshold = _std(dists) * factor

    for leaf in node.iter_leaves():
        if node.distance_to(leaf) > threshold:
            leaves_to_remove.add(leaf)

    node.remove_nodes(leaves_to_remove)
    return leaves_to_remove

def trim_zero_len_branches(node, min_len=1e-7):
    leaves_to_remove = set()

    dists = node.distances()
    for pair in dists:
        if dists[pair] < min_len:
            leaves_to_remove.add(pair[0], pair[1])

    return leaves_to_remove

def collapse_nodes(node, threshold):
    """
    Takes a TreeNode object and an integer or float support value as an input.
    Collapse all branches with a support value below the provided threshold into
    polytomies.
    """
    count = 0

    for branch in node.iter_branches():
        if not branch.support:
            continue

        if branch.support < threshold:
            branch.collapse()
            count += 1

    return count
