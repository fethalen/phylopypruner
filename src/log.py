# pylint: disable=too-many-instance-attributes
# pylint: disable=too-many-branches

"""
Store information about performed operations.
"""

from __future__ import print_function

class Log(object):
    """
    A record of a single run.
    """

    def __init__(self, version, msa, tree, settings):
        self._version = version
        self._msa = msa
        self._tree = tree
        self._msa_file = settings.fasta_file
        self._tree_file = settings.nw_file
        self._outgroup = settings.outgroup
        self._prune_paralogs = bool(settings.prune)
        self._sequences = len(list(self._tree.iter_leaves()))
        self._taxa = len(set(list(self._tree.iter_otus())))
        self._collapsed_nodes = 0
        self._trimmed_seqs = []
        self._lbs_removed = []
        self._monophylies_masked = []
        self._orthologs = []
        self._paralogs = []

    @property
    def version(self):
        """
        The version number used for this run.
        """
        return self._version

    @version.setter
    def version(self, value):
        self._version = value

    @property
    def msa(self):
        """
        The MultipleSequenceAlignment object used for this run.
        """
        return self._msa

    @msa.setter
    def msa(self, value):
        self._msa = value

    @property
    def tree(self):
        """
        The root of the TreeNode object used for this run.
        """
        return self._msa

    @tree.setter
    def tree(self, value):
        self._tree = value

    @property
    def msa_file(self):
        """
        The name of the MSA file used in this run.
        """
        return self._msa_file

    @msa_file.setter
    def msa_file(self, value):
        self._msa_file = value

    @property
    def tree_file(self):
        """
        The name of the Newick file used in this run.
        """
        return self._tree_file

    @tree_file.setter
    def tree_file(self, value):
        self._tree_file = value

    @property
    def outgroup(self):
        """
        A list of OTUs used as an outgroup in this run.
        """
        return self._outgroup

    @outgroup.setter
    def outgroup(self, value):
        self._outgroup = value

    @property
    def sequences(self):
        """
        Number of unique sequences used in this run.
        """
        return self._sequences

    @sequences.setter
    def sequences(self, value):
        self._sequences = value

    @property
    def trimmed_seqs(self):
        """
        A list of sequences that were deleted due to being to short.
        """
        return self._trimmed_seqs

    @property
    def prune_paralogs(self):
        """
        True if a prune paralog method was specified for this run.
        """
        return self._prune_paralogs

    @prune_paralogs.setter
    def prune_paralogs(self, value):
        self._prune_paralogs = value

    @trimmed_seqs.setter
    def trimmed_seqs(self, value):
        self._trimmed_seqs = value

    @property
    def monophylies_masked(self):
        """
        A list of sequences that were removed during monophyletic masking.
        """
        return self._monophylies_masked

    @monophylies_masked.setter
    def monophylies_masked(self, value):
        self._monophylies_masked = value

    @property
    def collapsed_nodes(self):
        """
        Number of collapsed nodes.
        """
        return self._collapsed_nodes

    @collapsed_nodes.setter
    def collapsed_nodes(self, value):
        self._collapsed_nodes = value

    @property
    def taxa(self):
        """
        Number of OTUs in the tree.
        """
        return self._taxa

    @taxa.setter
    def taxa(self, value):
        self._taxa = value

    @property
    def lbs_removed(self):
        """
        A list of long branches (LBs) that were removed during the lb-prune
        stage.
        """
        return self._lbs_removed

    @lbs_removed.setter
    def lbs_removed(self, value):
        self._lbs_removed = value

    @property
    def orthologs(self):
        "A list of TreeNode objects recovered as orthologs."
        return self._orthologs

    @orthologs.setter
    def orthologs(self, value):
        self._orthologs = value

    @property
    def paralogs(self):
        "A list of TreeNode objects recovered as paralogs."
        return self._paralogs

    @paralogs.setter
    def paralogs(self, value):
        self._paralogs = value

    def report(self, verbose):
        """
        Print a report of the records in this log.
        """
        print("MSA:\t\t\t\t\t{}".format(self.msa_file))
        print("tree:\t\t\t\t\t{}".format(self.tree_file))

        if self.outgroup:
            if len(self.outgroup) == 1:
                print("outgroup:\t\t\t\t{}".format(self.outgroup[0]))
            else:
                print("outgroups:\t\t\t\t", end="")
                for index, otu in enumerate(self.outgroup):
                    if index == len(self.outgroup) - 1:
                        print("{}".format(otu))
                    else:
                        print("{}, ".format(otu), end="")

        print("# of sequences:\t\t\t\t{}".format(self.sequences))
        print("# of OTUs:\t\t\t\t{}".format(self.taxa))
        print("# of short sequences removed:\t\t{}".format(len(self.trimmed_seqs)))
        print("# of long branched sequences removed\t{}".format(
            len(self.lbs_removed)))

        if self.monophylies_masked:
            print("# of monophylies masked:\t\t{}".format(
                len(self.monophylies_masked)))

        if self.collapsed_nodes:
            print("# of nodes collapsed into polytomies:\t{}".format(
                self.collapsed_nodes))

        if verbose:
            if self.trimmed_seqs:
                print("\nshort sequences removed:")
                for trimmed_seq in self.trimmed_seqs:
                    print("  {}".format(trimmed_seq))

            if self.lbs_removed:
                print("\nlong branched sequences removed:")
                for long_branch in self.lbs_removed:
                    print("  {}".format(long_branch))

        seen = set()
        if self.paralogs:
            print("\nOTUs with paralogs:")
            for paralog in self.paralogs:
                if not paralog.otu() in seen:
                    print("  {}".format(paralog.otu()))
                    seen.add(paralog.otu())

        if self.orthologs:
            for index, subtree in enumerate(self.orthologs):
                leaf_count = len(list(subtree.iter_leaves()))
                print("\northologous group #{}:\t\t\t\t".format(index + 1))
                print("  # of OTUs: \t\t\t\t{}".format(leaf_count))
                subtree.view()
        elif self.prune_paralogs:
            print("no orthologs were recovered")

    def to_txt(self, filename):
        """
        Takes a filename as an input and writes the records of this log to a
        plain text file.
        """
        pass

    def to_csv(self, filename):
        """
        Takes a filename as an input and writes the records in this log to a
        CSV file to the provided path.
        """
        pass

    def merge_logs(self, log):
        """
        Takes a Log object as an input and merges the fields within that log
        with this log.
        """
        pass
