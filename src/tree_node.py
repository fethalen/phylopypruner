"""
Module for working with a phylogenetic tree.
"""

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

    def traverse(self):
        """
        Returns an iterator object that includes all nodes related to this one.
        """
        nodes = self.children
        while nodes:
            node = nodes.pop()
            yield node
            for child in node.children:
                nodes.append(child)

    def iter_leaves(self):
        """
        Returns an iterator object that includes all leaves within this node.
        """
        for node in self.traverse():
            if node.is_leaf():
                yield node
