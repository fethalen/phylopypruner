"""
Methods for rooting a phylogenetic tree.
"""

import re

def outgroup(tree, outgroups):
    """
    Takes a TreeNode object and a list of outgroups as an input. If all
    provided are present within the tree and if they form a unique monophyletic
    group, then go ahead and root the tree using the ancestral node that
    includes all the outgroups.
    """
    if not tree.outgroups_present(outgroups):
        return

    if tree.repetetive_outgroups(outgroups):
        return

    if len(outgroups) == 1:
        outgroup_otu = outgroups[0]
        for leaf in tree.iter_leaves():
            otu = re.split(r"\||@", leaf.name)[0]
            if otu == outgroup_otu:
                rerooted_tree = tree.reroot(leaf)

    for branch in tree.iter_branches():
        if branch.is_monophyletic_outgroup(outgroups):
            for leaf in branch.iter_leaves():
                print(leaf.name)

    return rerooted_tree

def _longest_tip_to_tip_dist(node):
    """
    Takes a TreeNode object as an input. Calculates the tip to tip distances
    between the given node and the rest of the leaves within the tree and
    returns the longest distance.
    """

def midpoint(tree):
    """
    Takes a TreeNode object as an input, calculates the tip to tip distances
    and then roots the tree halfway between the two longest tips.
    """
    pass
    # longest_a, longest_b = None
    #
    # for node in tree.traverse_postorder():
    #     if not longest_a:
    #         longest_a = node.distance_to(

def molecular_clock(tree):
    """
    Takes a TreeNode object as an input ...
    """
    pass
