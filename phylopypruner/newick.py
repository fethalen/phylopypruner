#!/usr/bin/env python

"""
Class for reading and writing to or from a Newick file.
"""

from __future__ import absolute_import
import re
from collections import Counter

def _validate_newick(nw_str):
    """
    Takes a Newick string as an input and performs a bunch of tests to verify
    that the input is in the correct format.
    """
    if not nw_str.startswith("(") and nw_str.endswith(";"):
        raise AssertionError("File should start with \"(\" end with ;")
    count = Counter(nw_str)
    if count["("] != count[")"]:
        raise AssertionError("Unbalanced parentheses do not match.")
    return

def _read_label(label, node, interior_node):
    """
    Takes a node label, current working node and a boolean, signifying whether
    if the label belongs to an interior node or not, as an input. Adds the
    information to the current node, or a new node, depending on if it is an
    interior node or not.
    """
    label = label.split(")")[0]
    if ":" in label:
        label, branch_len = label.split(":")
    else:
        branch_len = 0 # No branch length is assigned, use default.

    if interior_node:
        node.dist = float(branch_len)
        if isinstance(label, str):
            node.name = label
        if isinstance(label, float) or isinstance(label, int):
            node.support = float(label)
    else:
        # Node label belongs to leaf, assign data to a new node.
        node.add_child(None, label, branch_len)
    return

def read(path, root=None):
    """
    Parse the Newick file contained in path and yield the root node of a
    TreeNode object.
    """
    with open(path, "r") as nw_file:
        nw_str = "".join(list(nw_file)) # merge lines to single string
        nw_str = re.sub("[\t\n\r ]", "", nw_str) # remove unnecessary whitespace
        _validate_newick(nw_str)

        # Create a new node if no node has been assigned.
        if not root:
            from phylopypruner.tree_node import TreeNode
            root = TreeNode()

        node = None
        for opening_paren in nw_str.split("(")[1:]:
            # Current node is root if no node has been assigned, else child
            # of previous node.
            node = root if not node else node.add_child()

            for subpart in opening_paren.split(","): # do for each comma
                if subpart:
                    _read_label(subpart, node, False)
                for closing_paren in subpart.split(")")[1:]:
                    closing_paren = closing_paren.rstrip(";")
                    _read_label(closing_paren, node, True)
                    node = node.parent # up one level
        return root
