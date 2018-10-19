# pylint: disable=too-many-instance-attributes
# pylint: disable=too-many-public-methods
# pylint: disable=too-many-branches

"Module for working with a phylogenetic tree."

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
        return re.split(r"\||@", self.name)[0]

    def is_root(self):
        "Returns True if this node lacks a parent."
        return not self.parent

    def is_leaf(self):
        "Returns True if this node has no children."
        return len(self.children) is 0

    def is_monophyletic(self):
        "Returns True if all leaves in this node belong to a single OTU."
        return len(set(list(self.iter_otus()))) is 1

    def is_polytomy(self):
        "Returns True if this node has two or more children."
        return len(self.children) > 2

    def is_bifurcating(self):
        "Returns True if this node has exactly two children."
        return len(self.children) is 2

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

        if not node in self.traverse_preorder():
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
        node.delete()
        root = TreeNode()
        root.add_child(node)
        root.add_child(sister)

        # remove unifurcations
        for node in root.traverse_preorder():
            parent = node.parent
            if not parent:
                continue

            if len(parent) is 1 and not parent.is_root():
                # Case where a child was removed from a bifurcating node; replace
                # the node with the child.
                child = parent.children[0]
                parent.remove_child(child)
                parent.name = child.name
                parent.dist = parent.dist + child.dist
                parent.support = child.support
                parent.children = child.children

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
        """
        Takes a set of strings, where each string is the name of a node to
        remove, as an input. Iterate over the nodes within this TreeNode object
        and remove a node if it matches a node within the set of nodes.
        """
        while nodes_to_remove:
            match = False

            for node in self.traverse_preorder():
                if node in nodes_to_remove:
                    if not node.is_root():
                        match = True
                        node.delete()
                        nodes_to_remove.remove(node)
                    break

            if not match:
                # no node to remove could be found
                break

        for leaf in self.iter_leaves():
            if leaf.name is None:
                leaf.delete()

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
        parent.remove_child(self)
        # self.parent = None

        if len(parent) is 1 and not parent.is_root():
            # Case where a child was removed from a bifurcating node; replace
            # the node with the child.
            child = parent.children[0]
            parent.remove_child(child)
            parent.name = child.name
            parent.dist = parent.dist + child.dist
            parent.support = child.support
            parent.children = child.children

        return self

    def outgroups_present(self, outgroups):
        """
        Takes a list of OTUs as an input. Returns true if all OTUs within
        outgroups are present within this TreeNode object.
        """
        for outgroup_otu in outgroups:
            found = False

            for otu in self.iter_otus():
                if otu == outgroup_otu:
                    found = True

            if not found:
                return False

        return True

    def repetetive_outgroups(self, outgroups):
        """
        Takes a list of OTUs as an input. Returns true if this TreeNode object
        contains multiple instances of the same OTU.
        """
        seen = set()

        for outgroup in self.iter_otus():
            if outgroup in seen:
                return True
            if outgroup in outgroups:
                seen.add(outgroup)

        return False

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
        """
        Takes a list of OTUs as an input. Returns true if this node only
        contains non-repetetive OTUs from the provided list.
        """
        if not isinstance(outgroups, list):
            raise TypeError("{} is not a list".format(outgroups))

        seen = set()

        for outgroup in self.iter_otus():
            if not outgroup in outgroups or outgroup in seen:
                return False
            seen.add(outgroup)

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

        self.delete()

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
        (lines, mid) = self._asciiArt(show_internal=show_internal,
                                      compact=compact, attributes=attributes)
        return '\n'+'\n'.join(lines)

    def empty_leaves(self):
        "Returns True if there are leaves without a name within this node."
        for leaf in self.iter_leaves():
            if leaf.is_root():
                return False
            if not leaf.name:
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
