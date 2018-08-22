# pylint: disable=too-many-instance-attributes
# pylint: disable=too-many-public-methods

"""
Module for working with a phylogenetic tree.
"""

import re
from collections import deque

class TreeNode(object):
    """
    Represents a node in a phylogenetic tree. Nodes may have one parent and one
    or more children.
    """

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
            print('child {} not found in node {}'.format(child, self))

    def remove_node(self, node):
        """
        Remove the provided TreeNode object from this node, if it exists. No
        operation is performed if this is the root node.
        """
        if not isinstance(node, TreeNode):
            raise TypeError("{} must be a TreeNode object".format(node))

        if node.is_root():
            raise AssertionError("{} can't be a root node".format(node))

        for node_item in self.traverse_preorder():
            if node_item is node:
                parent = node.parent

                if parent:
                    parent.remove_child(node)

                    if len(parent) is 0 and not parent.is_root():
                        remove_node(parent)
        return node

    def iter_sisters(self):
        """
        Return an iterator object that contains this node's sister nodes.
        """
        if self.parent:
            for child in self.parent.children:
                if child is not self:
                    yield child

    def is_root(self):
        """Returns True if this node lacks a parent."""
        return not self.parent

    def is_leaf(self):
        """Returns True if this node has no children."""
        return len(self.children) is 0

    def is_monophyletic(self):
        """Returns True if all leaves in this node belong to a single OTU."""
        return len(set(list(self.iter_otus()))) is 1

    def is_polytomy(self):
        """Returns True if this node has two or more children."""
        return len(self.children) > 2

    def is_bifurcating(self):
        "Returns True if this node has exactly two children."
        return len(self.children) is 2

    def reroot(self, node):
        """
        Takes a TreeNode object as an input. If the provided TreeNode object
        exists within this tree, then go ahead and reroot the tree using that
        node.
        """
        if not isinstance(node, TreeNode):
            raise TypeError("{} must be a Treeroot object".format(node))

        if not node in self.traverse_preorder():
            raise AssertionError("{} is not present within tree \
                    {}".format(node, self))

        root = None
        parent = node.parent

        index = 0
        while parent:
            index += 1

            if not root:
                root = TreeNode()
                working_node = root
                last = root.add_child(node)

            if not parent.is_root():
                working_node = working_node.add_child(child=None,
                                                      name=parent.name,
                                                      dist=parent.dist,
                                                      support=parent.support)

            for child in parent.children:
                if not child is last:
                    working_node.children.append(child)

            last = parent
            parent = parent.parent
            if parent:
                parent.remove_child(last)

        return root

    def outgroups_only(self, outgroups):
        """
        Takes a list of OTUs as an input. Returns true if this node only
        contains OTUs from the provided list.
        """
        if not isinstance(outgroups, list):
            raise TypeError("{} is not a list".format(outgroups))

        statement = True

        for outgroup in self.iter_otus():
            if not outgroup in outgroups:
                statement = False

        return statement

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

    def collapse(self):
        """
        Collapse this node into a polytomy; move all its children to its parent
        node and then remove the node itself.
        """
        parent = self.parent
        if not parent:
            raise AssertionError("{} does not have any parents".format(self))

        for child in self.children:
            parent.add_child(child)

        self.children = []
        remove_node(self)

    def distance_to(self, node):
        "Returns the distance between this node and another node."
        if not isinstance(node, TreeNode):
            raise TypeError("{} must be a TreeNode object".format(node))

        distance = 0

        if node is self:
            return distance

        distance += self.dist
        parent = self.parent

        # traverse this node
        while parent:
            if self and node in parent.traverse_preorder():
                # smallest node that includes both nodes under comparision
                smallest_subtree = parent
            if parent.dist:
                distance += parent.dist
            parent = parent.parent

        distance += node.dist
        parent = node.parent

        # traverse provided node
        while parent:
            if parent is smallest_subtree:
                return distance

            if parent.dist:
                distance += parent.dist

            parent = parent.parent

    def iter_leaves(self):
        """
        Returns an iterator object that includes all leaves within this node.
        """
        for node in self.traverse_preorder():
            if node.is_leaf():
                yield node

    def iter_branches(self):
        """
        Returns an iterator object that includes all branches within this node.
        """
        for node in self.traverse_preorder():
            if not node.is_leaf():
                yield node

    def iter_names(self):
        """
        Returns an iterator object that includes all names within this node.
        """
        for leaf in self.iter_leaves():
            yield leaf.name

    def iter_otus(self):
        """
        Returns an iterator object that includes all OTUs within this node.
        """
        for name in self.iter_names():
            otu = re.split(r"\||@", name)[0]
            yield otu

    def iter_identifiers(self):
        """
        Returns an iterator object that includes all identifiers in this node.
        """
        for name in self.iter_names():
            identifier = re.split(r"\||@", name)[1]
            yield identifier

    def view(self):
        "Print this tree in ASCII characters."
        print(self._get_ascii(compact=False, show_internal=False))

    def _asciiArt(self, char1='-', show_internal=True, compact=False, attributes=None):
        """
        Returns the ASCII representation of the tree.
        Code based on the PyCogent GPL project.
        """
        if not attributes:
            attributes = ["name"]
        node_name = ', '.join(map(str, [getattr(self, v) for v in attributes if hasattr(self, v)]))

        LEN = max(3, len(node_name) if not self.children or show_internal else 3)
        PAD = ' ' * LEN
        PA = ' ' * (LEN-1)
        if not self.is_leaf():
            mids = []
            result = []
            for c in self.children:
                if len(self.children) == 1:
                    char2 = '/'
                elif c is self.children[0]:
                    char2 = '/'
                elif c is self.children[-1]:
                    char2 = '\\'
                else:
                    char2 = '-'
                (clines, mid) = c._asciiArt(char2, show_internal, compact, attributes)
                mids.append(mid+len(result))
                result.extend(clines)
                if not compact:
                    result.append('')
            if not compact:
                result.pop()
            (lo, hi, end) = (mids[0], mids[-1], len(result))
            prefixes = [PAD] * (lo+1) + [PA+'|'] * (hi-lo-1) + [PAD] * (end-hi)
            mid = int((lo + hi) / 2)
            prefixes[mid] = char1 + '-'*(LEN-2) + prefixes[mid][-1]
            result = [p+l for (p,l) in zip(prefixes, result)]
            if show_internal:
                stem = result[mid]
                result[mid] = stem[0] + node_name + stem[len(node_name)+1:]
            return (result, mid)
        else:
            return ([char1 + '-' + node_name], 0)

    def _get_ascii(self, show_internal=True, compact=False, attributes=None):
        """
        Returns a string containing an ascii drawing of the tree.
        :argument show_internal: includes internal edge names.
        :argument compact: use exactly one line per tip.
        :param attributes: A list of node attributes to shown in the
            ASCII representation.
        """
        (lines, mid) = self._asciiArt(show_internal=show_internal,
                                      compact=compact, attributes=attributes)
        return '\n'+'\n'.join(lines)

def remove_node(node):
    """
    Remove the provided node from its node, if it exists. If a node's
    parent node has no children after removal, then remove the parent node
    too. If the parent node is the root node, don't remove it regardless.
    """
    if not isinstance(node, TreeNode):
        raise TypeError("{} must be a TreeNode object".format(node))

    if node.is_root():
        raise AssertionError("can't remove a root node")

    parent = node.parent
    parent.remove_child(node)

    # In case one node is removed from a bifurcating node, that node is no
    # longer needed.
    if len(parent) is 1 and not parent.is_root():
        child = parent.children[0]
        parent.name = child.name
        parent.dist = parent.dist + child.dist
        parent.support = child.support
        parent.children = child.children
