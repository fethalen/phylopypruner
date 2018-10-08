# pylint: disable=too-many-instance-attributes
# pylint: disable=too-many-branches

"""
Store information about performed operations.
"""

from __future__ import print_function
from __future__ import absolute_import
from datetime import datetime
import os.path
from phylopypruner.msa import MultipleSequenceAlignment
from phylopypruner.tree_node import TreeNode
from phylopypruner.settings import Settings

TIMESTAMP = datetime.now().strftime("%Y-%m-%d")
ORTHO_STATS_PATH = "/{}_ppp_ortho_stats.csv".format(TIMESTAMP)
LOG_PATH = "/{}_ppp_run.log".format(TIMESTAMP)

class Log(object):
    """
    A record of a single run.
    """

    def __init__(self, version, msa=MultipleSequenceAlignment, tree=TreeNode,
                 settings=Settings):
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
        self._msas_out = []
        self._homology_tree = None
        self._masked_tree = None
        self._masked_tree_str = None
        self._settings = settings

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
        "True if a prune paralog method was specified for this run."
        return self._prune_paralogs

    @prune_paralogs.setter
    def prune_paralogs(self, value):
        self._prune_paralogs = value

    @property
    def msas_out(self):
        """
        A list of multiple sequence alignments (MSAs) that were added to the
        output.
        """
        return self._msas_out

    @msas_out.setter
    def msas_out(self, value):
        self._msas_out = value

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

    @property
    def homology_tree(self):
        "The tree as it looked before any operations were performed."
        return self._homology_tree

    @homology_tree.setter
    def homology_tree(self, value):
        self._homology_tree = value

    @property
    def masked_tree(self):
        """
        The tree after monophyletic masking and before the paralogy pruning
        stage.
        """
        return self._masked_tree

    @masked_tree.setter
    def masked_tree(self, value):
        self._masked_tree = value

    @property
    def masked_tree_str(self):
        """
        The tree after monophyletic masking and before the paralogy pruning
        stage as a string.
        """
        return self._masked_tree_str

    @masked_tree_str.setter
    def masked_tree_str(self, value):
        self._masked_tree_str = value

    @property
    def settings(self):
        "The Settings object that stores the settings used in this log."
        return self._settings

    @settings.setter
    def settings(self, value):
        self._settings = value

    def outgroups_to_str(self):
        """
        Returns a string that contains the outgroups that were used in this
        run.
        """
        if self.outgroup:
            if len(self.outgroup) == 1:
                return "outgroup:\t\t\t\t{}".format(self.outgroup[0])
            else:
                outgroups = "outgroups:\t\t\t\t"
                for index, otu in enumerate(self.outgroup):
                    if index == len(self.outgroup) - 1:
                        outgroups += otu
                    else:
                        outgroups += "{}, ".format(otu)
                return outgroups

    def paralogs_to_str(self):
        "Returns a string that contains the paralogs found in this run."
        unique_paralogs = set()
        seen = set()
        paralog_str = "paralogous OTUs: "
        if self.paralogs:
            for paralog in self.paralogs:
                if not paralog.otu() in seen:
                    unique_paralogs.add(paralog.otu())
                    seen.add(paralog.otu())
            for index, paralog in enumerate(unique_paralogs):
                if index == len(unique_paralogs) - 1:
                    paralog_str += paralog
                else:
                    paralog_str += "{}, ".format(paralog)
            return paralog_str
        else:
            return paralog_str + " none"

    def msas_out_to_str(self):
        """
        Returns a string that contains the name of the files that were written
        for this run.
        """
        msa_out_str = str()
        for msa_out in self.msas_out:
            if msa_out is self.msas_out[0]:
                msa_out_str += "\n"
            if msa_out is self.msas_out[-1]:
                msa_out_str += "wrote: {}\n".format(str(msa_out))
            else:
                msa_out_str += "wrote: {}".format(str(msa_out))
        return msa_out_str

    def orthologs_to_str(self):
        """
        Returns a string that contains statistics for the orthologs found in
        this run.
        """
        ortho_str = ""
        if self.orthologs:
            for index, subtree in enumerate(self.orthologs):
                if subtree:
                    leaf_count = len(list(subtree.iter_leaves()))
                    ortho_str += "\northologous group #{}:\t\t\t\t".format(index + 1)
                    ortho_str += "\n  # of sequences:\t{}".format(leaf_count)
                    ortho_str += "\n{}".format(subtree.view())
        else:
            ortho_str = "no orthologs were recovered"
        return ortho_str

    def report(self, dir_out):
        "Print a report of the records in this log."
        report = """
MSA:\t\t\t\t\t{}
tree:\t\t\t\t\t{}
{}
# of sequences:\t\t\t\t{}
# of OTUs:\t\t\t\t{}
# of short sequences removed:\t\t{}
# of long branched sequences removed:\t{}
# of monophylies masked:\t\t{}
# of nodes collapsed into polytomies:\t{}
{}
\ninput tree:\n{}
\ntree before paralogy pruning:\n{}
{}
{}""".format(self.msa_file, self.tree_file, self.outgroups_to_str(),
             self.sequences, self.taxa, len(self._trimmed_seqs),
             len(self.lbs_removed), len(self.monophylies_masked),
             self.collapsed_nodes, self.paralogs_to_str(),
             self.homology_tree, self.masked_tree_str,
             self.orthologs_to_str(), self.msas_out_to_str())

        with open(dir_out + LOG_PATH, "a") as log_file:
            log_file.write(report)

        for ortholog in self.orthologs:
            self.to_csv(dir_out, ortholog)

    def to_csv(self, dir_out, ortholog):
        """
        Takes a filename as an input and writes the records in this log to a
        CSV file to the provided path.
        """
        with open(dir_out + ORTHO_STATS_PATH, "a") as stats_file:
            for ortholog in self.msas_out:
                row = "{};{};{};{};{};{};{}\n".format(
                    os.path.basename(str(ortholog)),
                    len(ortholog),
                    avg_seq_len(ortholog),
                    shortest_sequence(ortholog),
                    longest_sequence(ortholog),
                    missing_data(ortholog),
                    cat_alignment(ortholog))
                stats_file.write(row)

    def msa_out_path(self, dir_out, index=""):
        msa_in_path = str(self.msa)
        msa_in_filename = os.path.basename(msa_in_path)
        basename, extension = os.path.splitext(msa_in_filename)

        orthologs_dir = "{}/{}_orthologs".format(dir_out, TIMESTAMP)
        if not os.path.isdir(orthologs_dir):
            os.makedirs(orthologs_dir)

        if index:
            index = "_{}".format(index)

        return("{}/{}_pruned{}{}".format(orthologs_dir, basename, index, extension),
               extension)

    def get_msas_out(self, dir_out):
        """
        Takes an output directory as an input and generates a list of
        MultipleSequenceAlignments based on the TreeNode objects within
        this Log object's list of orthologs. The MSAs are stored in this logs
        list of MSAs called 'msas_out'.
        """
        for index, ortholog in enumerate(self.orthologs):
            if len(self.orthologs) is 1:
                msa_out_path, extension = self.msa_out_path(dir_out)
            else:
                msa_out_path, extension = self.msa_out_path(dir_out, index)

            msa_out = MultipleSequenceAlignment(msa_out_path, extension)

            for leaf in ortholog.iter_leaves():
                try:
                    msa_out.add_sequence(self.msa.get_sequence(leaf.name))
                except:
                    pass

            if len(msa_out) > 0:
                self.msas_out.append(msa_out)

def avg_seq_len(msa):
    """
    Takes a MultipleSequenceAlignment object as an input and returns the
    average sequence length of the sequences within it.
    """
    seq_lens = 0
    sequences = 0
    for sequence in msa.sequences:
        sequences += 1
        seq_lens += len(sequence.ungapped())
    if sequences > 0:
        return round(seq_lens / sequences, 1)
    else:
        return 0

def missing_data(msa):
    "Returns the percent missing data within the ortholog."
    sequences = 0
    pct_missing = 0.0
    for sequence in msa.sequences:
        sequences += 1
        ungapped = len(sequence.ungapped())
        missing = len(sequence) - ungapped
        if len(sequence) > 0:
            pct_missing += float(missing) / float(len(sequence))
    if sequences > 0:
        return round(pct_missing / sequences, 3) * 100
    else:
        return 0

def cat_alignment(msa):
    "Returns the length of the alignment."
    seq_lens = 0
    # It doesn't matter which sequence we picked since they're
    # aligned and missing positions are denoted by gaps.
    sequence = msa.sequences[0]
    seq_lens += len(sequence)
    return seq_lens

def shortest_sequence(msa):
    """
    Returns the shortest sequence in the provided MultipleSequenceAlignment
    object.
    """
    shortest = None
    for sequence in msa.sequences:
        if not shortest or shortest > len(sequence.ungapped()):
            shortest = len(sequence.ungapped())
    return shortest

def longest_sequence(msa):
    """
    Returns the longest sequence in the provided MultipleSequenceAlignment
    object.
    """
    longest = None
    for sequence in msa.sequences:
        if not longest or longest < len(sequence.ungapped()):
            longest = len(sequence.ungapped())
    return longest

def _file_out(path, dir_out=None, index=""):
    """
    Takes the path to an MSA and an optional path to a directory as an input.
    Extracts the base name and extension from the provided path. If no directory
    has been specified, then the directory is also extracted from the path.
    Returns the path to a file in the following format:
      <dir_out>/<basename>_pruned<extension>

    If an index has been provided, then the output will be in the following
    format:
      <dir_out>/<basename>_pruned_<index><extension>
    """
    if not directory:
        directory = os.path.dirname(path)
    filename = os.path.basename(path)
    basename, extension = os.path.splitext(filename)

    dir_out = directory + "/{}_orthologs".format(TIMESTAMP)
    if not os.path.isdir(dir_out):
        os.makedirs(dir_out)

    if index:
        index = "_{}".format(index)

    return "{}/{}_pruned{}{}".format(dir_out, basename, index, extension)
