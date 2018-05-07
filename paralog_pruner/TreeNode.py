#!/usr/bin/env python

""" Class for representing a node in a phylogenetic tree.
"""

class TreeNode(object):
    """ Represents a node that can be used to store a phylogenetic tree.
    """
    def __init__(self, child=None, name='', dist=0.0, support=1.0):
        self._name = name
        self._children = []
        self._parent = None
        self._dist = dist
        self._support = support

        self.name = name

    def __str__(self):
        return self.name

    def __len__(self):
        """ Length of node is its number of children. """
        return len(self.children)

    def __nonzero__(self):
        return True

    def __bool__(self):
        """ Python3 __nonzero__ equivalent. """
        return True

    def _get_children(self):
        return self._children

    def _set_children(self, value):
        if isinstance(value, list):
            self._children = value

    def _get_parent(self):
        return self._parent

    def _set_parent(self, value):
        self._parent = value

    def _get_dist(self):
        return self._dist

    def _set_dist(self, value):
        self._dist = float(value)

    def _get_support(self):
        return self._support

    def _set_support(self, value):
        self._support = float(value)

    def add_child(self, child=None, name='', dist=0.0, support=1.0):
        """ Adds a new child to this node.
        """
        if not child:
            child = self.__class__()

        child.name = name
        child.dist = dist
        child.support = support

        self.children.append(child)
        child.parent = self
        return child

    def remove_child(self, child):
        """ Takes the name of a child in this node and removes it, if it
        exists.
        """
        try:
            self._children.remove(child)
        except ValueError:
            print('Child not found.')

    dist = property(fget=_get_dist, fset=_set_dist)
    support = property(fget=_get_support, fset=_set_support)
    parent = property(fget=_get_parent, fset=_set_parent)
    children = property(fget=_get_children, fset=_set_children)

    def is_root(self):
        """ Returns true if this node lacks a parent. """
        return not self.parent

    def is_leaf(self):
        """ Returns True if this node has no children. """
        return len(self.children) is 0

    def traverse(self):
        """ Returns an iterator object that includes all nodes related to this
        one.
        """
        nodes = self.children
        while nodes:
            node = nodes.pop()
            yield node
            for child in node.children:
                nodes.append(child)

    def iter_leaves(self):
        """ Returns an iterator object that includes all leaves within this
        node.
        """
        for node in self.traverse():
            if node.is_leaf():
                yield node
