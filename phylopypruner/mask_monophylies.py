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

    return removed

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
        if branch.is_monophyletic():
            # branch is monophyletic, find the longest sequence
            for leaf in branch.iter_leaves():
                seq_len = _sequence_len(msa, leaf.name)

                if seq_len > longest_seq:
                    longest_seq = seq_len
                    keep = leaf.name

            removed += _mask(branch, keep)
        else:
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

        longest_seq = 0
        keep = None

    for branch in node.iter_branches():
        if branch.parent:
            monophyletic_sister = _monophyletic_sister(branch)
            if monophyletic_sister:
                removed.append(monophyletic_sister)

    new_root = _get_new_root(node)
    return new_root, removed

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

def _monophyletic_branch(node):
    keep = _smallest_distance(node)
    return node.leaves_except(keep)

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

    shortest = None
    leaves_to_remove = set()

    for otu in multiples:
        for leaf in multiples[otu]:
            distance = leaf.distance_to(node)
            if not shortest or distance < shortest:
                keep = leaf.name
                shortest = distance

        for leaf in multiples[otu]:
            if leaf.name is not keep:
                leaves_to_remove.add(leaf)
    return leaves_to_remove

def pairwise_distance(node):
    leaves_to_remove = set()
    masked = set()
    monophylies = False

    for branch in node.iter_branches():
        monophyletic_sister = _monophyletic_sister(branch)
        if monophyletic_sister:
            leaves_to_remove.add(monophyletic_sister)
        elif branch.is_monophyletic():
            leaves_to_remove.update(_monophyletic_branch(branch))
        else:
            # one or more repetetive OTUs within the polytomy
            leaves_to_remove.update(_monophyletic_polytomy(branch))
    monophylies = bool(leaves_to_remove)

    # new monophylies may arise when we get rid of others
    while monophylies:
        masked.update(leaves_to_remove)
        node.remove_nodes(leaves_to_remove)
        leaves_to_remove = set()

        for branch in node.iter_branches():
            monophyletic_sister = _monophyletic_sister(branch)
            if monophyletic_sister:
                leaves_to_remove.add(monophyletic_sister)
            elif branch.is_monophyletic():
                leaves_to_remove.update(_monophyletic_branch(branch))
            else:
                # one or more repetetive OTUs within the polytomy
                leaves_to_remove.update(_monophyletic_polytomy(branch))
        monophylies = bool(leaves_to_remove)

    return node, masked

def pairwise_distance(node):
    print("before monophyly masking: {}".format(len(set(node.iter_leaves()))))
    masked = set()
    leaves_to_remove = set()
    for branch in node.iter_branches():
        if branch.is_monophyletic() and len(set(branch.iter_leaves())) is 2:
            leaves_to_remove = _monophyletic_branch(branch)
            break
    monophylies = bool(leaves_to_remove)

    while monophylies:
        masked.update(leaves_to_remove)
        for leaf in leaves_to_remove:
            leaf.delete()
        # node.remove_nodes(leaves_to_remove)
        leaves_to_remove = set()
        for branch in node.iter_branches():
            no_of_leaves = len(set(branch.iter_leaves()))

            if branch.is_monophyletic():
                if no_of_leaves is 2:
                    # repetetive OTUs within a bifurcating node with two leaves
                    leaves_to_remove = _monophyletic_branch(branch)
                elif no_of_leaves > 2:
                    leaves_only = True
                    for child in branch.children:
                        if not child.is_leaf():
                            leaves_only = False
                    if leaves_only:
                        # branch is a polytomy, where each child is a leaf and
                        # belongs to the same OTU
                        leaves_to_remove = _monophyletic_branch(branch)
            else:
                # look for one or more repetetive OTUs within the polytomy
                leaves_to_remove.update(_monophyletic_polytomy(branch))
            if leaves_to_remove:
                break
        monophylies = bool(leaves_to_remove)

    print("sequences masked: {}".format(len(masked)))
    print("after monophyly masking: {}".format(len(set(node.iter_leaves()))))

    # for branch in node.iter_branches():
    #     child = _monophyletic_sister(branch)
    #     if child:
    #         child.parent.remove_child(child)
    #         masked.add(child)
    node = _get_new_root(node)
    return node, masked
