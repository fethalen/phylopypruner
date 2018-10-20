"Module for storing and retrieving a set of options."

class Settings(object):
    """
    A list of settings used in a single run.
    """
    def __init__(self, arguments=None):
        self._fasta_file = arguments.msa
        self._nw_file = arguments.tree
        self._min_taxa = arguments.min_taxa
        self._min_len = arguments.min_len
        self._min_support = arguments.min_support
        self._trim_lb = arguments.trim_lb
        self._outgroup = arguments.outgroup
        self._exclude = arguments.exclude
        self._root = arguments.root
        self._mask = arguments.mask
        self._prune = arguments.prune
        self._trim_divergent = arguments.trim_divergent

    @property
    def fasta_file(self):
        "Path to the FASTA file."
        return self._fasta_file

    @fasta_file.setter
    def fasta_file(self, value):
        self._fasta_file = value

    @property
    def nw_file(self):
        "Path to the Newick tree file."
        return self._nw_file

    @nw_file.setter
    def nw_file(self, value):
        self._nw_file = value

    @property
    def min_taxa(self):
        "Minimum number of OTUs allowed in the output."
        return self._min_taxa

    @min_taxa.setter
    def min_taxa(self, value):
        self._min_taxa = value

    @property
    def min_len(self):
        "Minimum number of positions allowed in a sequence."
        return self._min_len

    @min_len.setter
    def min_len(self, value):
        self._min_len = value

    @property
    def min_support(self):
        "Minimum support value allowed in a tree."
        return self._min_support

    @min_support.setter
    def min_support(self, value):
        self._min_support = value

    @property
    def trim_lb(self):
        "Integer factor used for getting rid of long branches."
        return self._trim_lb

    @trim_lb.setter
    def trim_lb(self, value):
        self._trim_lb = value

    @property
    def outgroup(self):
        "List of outgroups."
        return self._outgroup

    @outgroup.setter
    def outgroup(self, value):
        self._outgroup = value

    @property
    def root(self):
        "Rooting method used."
        return self._root

    @root.setter
    def root(self, value):
        self._root = value

    @property
    def mask(self):
        "Masking method used."
        return self._mask

    @mask.setter
    def mask(self, value):
        self._mask = value

    @property
    def prune(self):
        "Pruning method used."
        return self._prune

    @prune.setter
    def prune(self, value):
        self._prune = value

    @property
    def exclude(self):
        "A list of taxa to exclude."
        return self._exclude

    @exclude.setter
    def exclude(self, value):
        self._exclude = value

    @property
    def trim_divergent(self):
        "Trim divergent sequence threshold."
        return self._trim_divergent

    @trim_divergent.setter
    def trim_divergent(self, value):
        self._trim_divergent = value
