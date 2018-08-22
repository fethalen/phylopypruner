"""
Tools for filtering tree nodes and sequences.
"""

from tree_node import remove_node

def _is_short_sequence(sequence, treshold):
    """
    Return true if the provided sequence is shorter than the provided treshold.
    """
    return len(sequence.ungapped()) < treshold

def too_few_otus(tree, treshold):
    """
    Return true if the provided tree node object have fewer OTUs than the
    provided treshold.
    """
    if len(list(tree.iter_otus())) < treshold:
        return True

def trim_short_seqs(msa, tree, treshold):
    """
    Takes a TreeNode object, an MSA object and a treshold as an input. Remove
    sequences that are shorter than the provided treshold from both the MSA and
    the tree.
    """
    trimmed_seqs = []

    for leaf in tree.iter_leaves():
        match = msa.get_sequence(leaf.name)

        if match:
            sequence = match

            if _is_short_sequence(sequence, treshold):
                remove_node(leaf)
                trimmed_seqs.append(leaf)

    # TODO: nothing is appended!
    return trimmed_seqs

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
    """ Return the population standard deviation of data.
    """
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
    for leaf in node.iter_leaves():
        dists.append(leaf.dist)
    treshold = _std(dists) * factor

    for leaf in node.iter_leaves():
        if leaf.dist > treshold:
            remove_node(leaf)
            yield leaf

def collapse_nodes(node, treshold):
    """
    Takes a TreeNode object and an integer or float support value as an input.
    Collapse all branches with a support value below the provided treshold into
    polytomies.
    """
    count = 0

    for branch in node.iter_branches():
        if not branch.support:
            continue

        if branch.support < treshold:
            branch.collapse()
            count += 1

    return count
