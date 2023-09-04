"Module for working with a phylogenetic tree."

from __future__ import absolute_import
import re
from collections import deque
from . import report

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
        "Length of node is its number of children."
        return len(self.children)

    def __nonzero__(self):
        return True

    def __bool__(self):
        "Python3 __nonzero__ equivalent."
        return True

    @property
    def name(self):
        "Return the name of this node."
        return self._name

    @name.setter
    def name(self, value):
        self._name = str(value)

    @property
    def children(self):
        "Return a list of the node's children."
        return self._children

    @children.setter
    def children(self, value):
        if isinstance(value, list):
            self._children = value

    @property
    def parent(self):
        "Return the node's parent, if any, or return None."
        return self._parent

    @parent.setter
    def parent(self, value):
        self._parent = value

    @property
    def dist(self):
        "This node's distance."
        return self._dist

    @dist.setter
    def dist(self, value):
        self._dist = float(value)

    @property
    def support(self):
        "The support value for this node. Default is 1.0."
        return self._support

    @support.setter
    def support(self, value):
        support_value = float(value)
        if support_value > 1 and support_value <= 100:
            support_value = support_value / 100
        self._support = support_value

    def add_child(self, child=None, name='', dist=0.0, support=1.0):
        "Adds a new child to this node."
        if not child:
            child = self.__class__()
            child.name = name
            child.dist = dist
            child.support = support

        self.children.append(child)
        child.parent = self
        return child

    def remove_child(self, child):
        "Remove the provided child, if it exists."
        try:
            self._children.remove(child)
        except ValueError:
            print('child {} not found in node {}'.format(child, self))

    def iter_sisters(self):
        "Return an iterator object that contains this node's sister nodes."
        if self.parent:
            for child in self.parent.children:
                if child is not self:
                    yield child

    def otu(self):
        "Returns the OTU to which this node belongs."
        return re.split(r"\||@|_", self.name)[0]

    def is_root(self):
        "Returns True if this node lacks a parent."
        return not self.parent

    def is_leaf(self):
        "Returns True if this node has no children."
        return len(self.children) == 0

    def is_monophyletic(self):
        "Returns True if all leaves in this node belong to a single OTU."
        return len(set(list(self.iter_otus()))) == 1

    def is_polytomy(self):
        "Returns True if this node has two or more children."
        return len(self.children) > 2

    def is_bifurcating(self):
        "Returns True if this node has exactly two children."
        return len(self.children) == 2

    def reroot(self, node):
        """Root this node at the provided node.

        Parameters
        ----------
        node : TreeNode object
            The TreeNode object that becomes the new root, must exist within
            this TreeNode object.

        Returns
        -------
        root : TreeNode object
            This tree rooted at the provided node.
        """
        if not isinstance(node, TreeNode):
            raise TypeError("{} must be a Treeroot object".format(node))

        if node not in self.traverse_preorder():
            raise AssertionError("{} is not present within tree \
                    {}".format(node, self))

        parent = node.parent
        ancestors = list()

        while parent:
            ancestors.append(parent)
            parent = parent.parent

        for index in range(len(ancestors) - 1):
            child = ancestors[index]
            parent = ancestors[index + 1]
            if child in parent.children:
                parent.remove_child(child)
                child.add_child(parent)

        sister = node.parent
        # node.delete()
        if node.parent:
            node.parent.remove_child(node)
        root = TreeNode()
        root.add_child(node)
        root.add_child(sister)

        # Remove monofurcating nodes.
        while self.has_monofurcations():
            self.prune_monofurcations()

        # remove empty leaves
        while root.empty_leaves():
            for leaf in root.iter_leaves():
                if not leaf.name:
                    leaf.delete()
                    break

        # workaround to make sure that child's parent is the same node as to
        # which the child belongs to
        for node in root.traverse_preorder():
            if not node.children:
                continue

            for child in node.children:
                child.parent = node

        self = root
        return self

    def remove_nodes(self, nodes_to_remove):
        """Takes a set of strings, where each string is the name of a node to
        remove. Iterate over the nodes within this TreeNode object and remove
        the node, if and only if its name matches one of the names within the
        provided set.

        Parameters
        ----------
        nodes_to_remove : set
            A set of strings where each string the name of a node that should
            be removed.

        Returns
        -------
        self : TreeNode
            This TreeNode object with nodes with names that matches the
            provided names removed
        """
        match = False
        while nodes_to_remove:
            match = False
            for node in self.traverse_preorder():
                if node in nodes_to_remove and not node.is_root():
                    match = True
                    try:
                        node.delete()
                    except:
                        report.error("can't delete node {}".format(node))
                    nodes_to_remove.remove(node)

            if not match:
                break

        # Remove monofurcating nodes due to deleted nodes.
        while self.has_monofurcations():
            self.prune_monofurcations()

        # remove empty leaves
        while self.empty_leaves() and not len(self) <= 1:
            for leaf in self.iter_leaves():
                if not leaf.name:
                    try:
                        leaf.delete()
                    except:
                        report.error("can't delete empty leaf")
                    break

        # workaround to make sure that child's parent is the same node as to
        # which the child belongs to
        for node in self.traverse_preorder():
            if not self.children:
                continue

            for child in node.children:
                child.parent = node

        return self

    def delete(self):
        """
        Remove this node from its parent, if it exists. If the parent is a
        bifurcating node that is left with a single child upon removal (i.e.,
        the parent is unifurcating upon removal), then replace the node with
        its child.
        """
        if self.is_root():
            raise AssertionError("can't remove a root node")

        parent = self.parent

        if not parent:
            return self

        if len(self.children) == 1:
            self.children[0].dist += self.dist
        elif len(self.children) > 1:
            parent.dist += self.dist

        for child in self.children:
            parent.add_child(child)

        parent.remove_child(self)

        if parent.is_monofurcating():
            parent.delete()

        return self

    def is_monofurcating(self):
        """Returns True if the number of children that belongs to this node is
        equal to 1 and this is not a root node.

        Returns
        -------
        True if the number of children in this node == 1, else False.
        """
        return not self.is_leaf() and\
            not self.is_root() and\
            len(self.children) < 2

    def has_monofurcations(self):
        """Returns True if there are any non-root nodes which are
        monofurcating, meaning that they have only one child.

        Returns
        -------
        True if any non-root node within this node is monofurcating, else
        False.
        """
        for branch in self.iter_branches():
            if branch.is_monofurcating():
                return True
        return False

    def prune_monofurcation(self):
        """Remove a monofurcating node, meaning a non-root node with a single
        child.

        Returns
        -------
        self : TreeNode object
            The monofurcating node, replaced with its child.
        """
        child = self.children[0]
        child.delete()
        return self

    def prune_monofurcations(self):
        """Remove multiple monofurcating nodes (a non-root with a single child)
        within this node.

        Returns
        -------
        self : TreeNode object
            This node without monofurcations.
        """
        for node in self.traverse_preorder():
            parent = node.parent
            if not parent:
                continue

            if parent.is_monofurcating():
                parent.prune_monofurcation()

        return self

    def outgroups_present(self, outgroups):
        """Takes a list of OTUs as an input and returns the subset of those
        OTUs that are present within this tree.

        Parameters
        ----------
        outgroups : list
            A list of outgroups to look for.

        Returns
        -------
        outgroups_in_node : list
            A list of the outgroups that are present within this TreeNode
            object.
        """
        outgroups_in_node = []

        for otu in self.iter_otus():
            if otu in outgroups:
                outgroups_in_node.append(otu)

        return outgroups_in_node

    def distances(self):
        """
        Return the distance for all pairs of leaves within this TreeNode
        object.
        """
        distances = dict()

        for leaf_a in self.iter_leaves():
            for leaf_b in self.iter_leaves():
                if leaf_a is leaf_b:
                    continue

                leaves_in_key = False
                for key in distances:
                    if leaf_a and leaf_b in key:
                        leaves_in_key = True

                if not leaves_in_key:
                    distances[(leaf_a, leaf_b)] = leaf_a.distance_to(leaf_b)

        return distances

    def is_monophyletic_outgroup(self, outgroups):
        """Takes a list of OTUs as an input. Returns true if this node only
        contains non-repetetive OTUs from the provided list.
        """
        if not isinstance(outgroups, list):
            raise TypeError("{} is not a list".format(outgroups))

        otus = set(self.iter_otus())

        # None of the outgroups is present within this TreeNode object.
        if not otus.intersection(outgroups):
            return False

        # Repetetive OTUs found within this TreeNode object.
        if len(set(outgroups)) < len(outgroups):
            return False

        for leaf in self.iter_leaves():
            if not leaf.otu() in outgroups:
                return False

        return True

    def traverse_preorder(self):
        "Traverse the tree, from the root to the leaves."
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

        parent.remove_child(self)
        # self.delete()

        return self

    def distance_to(self, node):
        "Returns the distance between this node and another node."
        if not isinstance(node, TreeNode):
            raise TypeError("{} must be a TreeNode object".format(node))

        distance = 0.0

        if node is self:
            return distance

        # We need a different approach if this node is a root or if the
        # provided node is a root.
        if self.is_root():
            root = self
            leaf = node
        elif node.is_root():
            root = node
            leaf = self
        else:
            root = None

        if root:
            distance += root.dist
            parent = leaf.parent

            while parent:
                distance += parent.dist
                parent = parent.parent
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

        while parent:
            if parent is smallest_subtree:
                return distance

            if parent.dist:
                distance += parent.dist

            parent = parent.parent

    def iter_leaves(self):
        "Returns an iterator object that includes all leaves within this node."
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
        "Returns an iterator object that includes all names within this node."
        for leaf in self.iter_leaves():
            yield leaf.name

    def iter_otus(self):
        "Returns an iterator object that includes all OTUs within this node."
        for name in self.iter_names():
            otu = re.split(r"\||@|_", name)[0]
            yield otu

    def iter_identifiers(self):
        """
        Returns an iterator object that includes all identifiers in this node.
        """
        for name in self.iter_names():
            identifier = re.search(r"[|@_]([^ ]*)", name).group(1)
            yield identifier

    def view(self):
        "Print this tree in ASCII characters."
        return self._get_ascii(compact=False, show_internal=False)

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
        (lines, _) = self._asciiArt(show_internal=show_internal,
                                      compact=compact, attributes=attributes)
        return '\n'+'\n'.join(lines)

    def empty_leaves(self):
        "Returns True if there are leaves without a name within this node."
        for leaf in self.iter_leaves():
            if not leaf.name and not leaf.is_root():
                return True
        return False

    def leaves_except(self, node):
        """
        Takes a TreeNode object as an input. Returns a set of leaves within
        this node that do not match the provided node.
        """
        leaves = set()
        for leaf in self.iter_leaves():
            if leaf is not node:
                leaves.add(leaf)
        return leaves

    def paralogs(self):
        """
        Returns the list of TreeNode objects within this node with an OTU that
        is present more than once within this node.
        """
        seen = set()
        flag = set()
        paralogs = list()

        for leaf in self.iter_leaves():
            otu = leaf.otu()
            if otu in seen:
                flag.add(otu)
            seen.add(otu)

        for leaf in self.iter_leaves():
            otu = leaf.otu()
            if otu in flag:
                paralogs.append(leaf)

        return paralogs
