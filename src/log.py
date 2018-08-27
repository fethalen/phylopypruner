# pylint: disable=too-many-instance-attributes
# pylint: disable=too-many-branches

"""
Store information about performed operations.
"""

import datetime

class Log(object):
    """
    A record of a single run.
    """

    def __init__(self, version, msa, tree, arguments):
        self._version = version
        self._msa = msa
        self._tree = tree
        self._msa_file = msa.filename
        self._tree_file = arguments.tree
        self._outgroup = arguments.outgroup
        self._sequences = len(list(self._tree.iter_leaves()))
        self._taxa = len(set(list(self._tree.iter_otus())))
        self._collapsed_nodes = 0
        self._trimmed_seqs = []
        self._lbs_removed = []
        self._monophylies_masked = []
        self._pruned_sequences = []

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
        return self._msa_file

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
    def pruned_sequences(self):
        """
        A list of sequences that were removed during the paralog-pruning stage.
        """
        return self._pruned_sequences

    @pruned_sequences.setter
    def pruned_sequences(self, value):
        self._pruned_sequences = value

    def report(self):
        """
        Print a report of the records in this log.
        """
        print(datetime.datetime.now().strftime("%A, %d. %B %Y %I:%M%p"))
        print("version number:\t\t\t\t{}".format(self.version))
        print("MSA:\t\t\t\t\t{}".format(self.msa_file))
        print("tree:\t\t\t\t\t{}".format(self.tree_file))

        if self.outgroup:
            if len(self.outgroup) == 1:
                print("outgroup:\t\t\t\t{}".format(self.outgroup[0]))
            else:
                print("outgroups: ")
                for otu in self.outgroup:
                    print("  {}".format(otu))

        print("# of sequences:\t\t\t\t{}".format(self.sequences))
        print("# of OTUs:\t\t\t\t{}".format(self.taxa))
        print("# of short sequences removed:\t\t{}".format(len(self.trimmed_seqs)))
        print("# of long branched sequences removed\t{}".format(
            len(self.lbs_removed)))

        if self.monophylies_masked:
            print("# of monophylies masked:\t\t{}".format(
                len(self.monophylies_masked)))
            print("# of nodes collapsed into polytomies:\t{}".format(
                self.collapsed_nodes))

        if self.trimmed_seqs:
            print("\nshort sequences removed:")
            for trimmed_seq in self.trimmed_seqs:
                print("  {}".format(trimmed_seq))

        if self.lbs_removed:
            print("\nlong branched sequences removed:")
            for long_branch in self.lbs_removed:
                print("  {}".format(long_branch))

        if self.pruned_sequences:
            for index, subtree in enumerate(self.pruned_sequences):
                leaf_count = len(list(subtree.iter_leaves()))
                print("\northolog # {}:\t\t\t\t".format(index + 1))
                print("  # of sequences in ortholog: \t\t{}".format(leaf_count))
                subtree.view()

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
