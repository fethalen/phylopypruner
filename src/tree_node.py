"""
Module for working with a phylogenetic tree.
"""

from collections import deque
import re

class TreeNode(object):
    """
    Represents a node in a phylogenetic tree. Nodes may have one parent and one
    or more children.
    """

    # pylint: disable=too-many-instance-attributes

    def __init__(self, name='', dist=0.0, support=1.0):
        self._name = name
        self._children = []
        self._parent = None
        self._dist = dist
        self._support = support

    def __str__(self):
        return self.name

    def __len__(self):
        """Length of node is its number of children."""
        return len(self.children)

    def __nonzero__(self):
        return True

    def __bool__(self):
        """Python3 __nonzero__ equivalent."""
        return True

    @property
    def name(self):
        """Return the name of this node."""
        return self._name

    @name.setter
    def name(self, value):
        self._name = str(value)

    @property
    def children(self):
        """Return a list of the node's children."""
        return self._children

    @children.setter
    def children(self, value):
        if isinstance(value, list):
            self._children = value

    @property
    def parent(self):
        """Return the node's parent, if any, or return None."""
        return self._parent

    @parent.setter
    def parent(self, value):
        self._parent = value

    @property
    def dist(self):
        """This node's distance (that is, it's branch length)."""
        return self._dist

    @dist.setter
    def dist(self, value):
        self._dist = float(value)

    @property
    def support(self):
        """The support value for this node. Default is 1.0."""
        return self._support

    @support.setter
    def support(self, value):
        self._support = float(value)

    def add_child(self, child=None, name='', dist=0.0, support=1.0):
        """Adds a new child to this node."""
        if not child:
            child = self.__class__()

        child.name = name
        child.dist = dist
        child.support = support

        self.children.append(child)
        child.parent = self
        return child

    def remove_child(self, child):
        """Remove the provided child, if it exists."""
        try:
            self._children.remove(child)
        except ValueError:
            print('Child not found.')

    def is_root(self):
        """Returns true if this node lacks a parent."""
        return not self.parent

    def is_leaf(self):
        """Returns True if this node has no children."""
        return len(self.children) is 0

    def traverse_preorder(self):
        """
        Traverse the tree, from the root to the leaves.
        """
        structure = deque()
        structure.append(self)
        while structure:
            node = structure.pop()
            yield node
            structure.extend(node.children)

    def traverse_inorder(self):
        pass

    def traverse_postorder(self):
        """
        Traverse the tree, starting from the leaves and working inward to the
        root.
        """
        structure_1, structure_2 = deque()
        structure_1.append(self)
        while structure_1:
            node = structure_1.pop()
            structure_2.append(node)
            structure_1.extend(node.children)
        while structure_2:
            yield structure_2.pop()

    def iter_leaves(self):
        """
        Returns an iterator object that includes all leaves within this node.
        """
        for node in self.traverse_preorder():
            if node.is_leaf():
                yield node

    def leaves(self):
        """
        Returns a list containing all leaves within this node.
        """
        return list(self.iter_leaves())

    def iter_otus(self):
        """
        Returns an iterator object that includes all OTUs within this node.
        OTUs are separated from a unique sequence identifier by a delimiter
        (either '@' or '|').
        """
        for leaf in self.iter_leaves():
            otu = re.split(r"\||@", leaf.name)[0]
            yield otu

    def otus(self):
        """
        Returns a list containing all OTUs within this node. OTUs are separated
        from a unique sequence identifier by a delimiter (either '@' or '|').
        """
        return list(self.iter_otus())
