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
        self._trim_zero_len = arguments.min_pdist
        self._outgroup = arguments.outgroup
        self._include = arguments.include
        self._exclude = arguments.exclude
        self._root = arguments.root
        self._mask = arguments.mask
        self._prune = arguments.prune
        self._trim_freq_paralogs = arguments.trim_freq_paralogs
        self._trim_divergent = arguments.trim_divergent
        self._jackknife = arguments.jackknife

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
    def include(self):
        "A list of taxa to include, despite being deemed problematic."
        return self._exclude

    @include.setter
    def include(self, value):
        self._include = value

    @property
    def trim_freq_paralogs(self):
        "Trim frequent paralogs factor."
        return self._trim_freq_paralogs

    @trim_freq_paralogs.setter
    def trim_freq_paralogs(self, value):
        self._trim_freq_paralogs = value

    @property
    def trim_divergent(self):
        "Trim divergent sequence threshold."
        return self._trim_divergent

    @trim_divergent.setter
    def trim_divergent(self, value):
        self._trim_divergent = value

    @property
    def jackknife(self):
        "True if jackknifing is set."
        return self._jackknife

    @jackknife.setter
    def jackknife(self, value):
        self._jackknife = value

    @property
    def trim_zero_len(self):
        "True if jackknifing is set."
        return self._trim_zero_len

    @trim_zero_len.setter
    def trim_zero_len(self, value):
        self._trim_zero_len = value

    def report(self, directory, log_path):
        "Generate a report of what settings were used."
        if directory:
            input_files = "Directory:\t{}".format(directory)
        else:
            input_files = "MSA:\t{}\n  Tree:\t{}".format(
                self.fasta_file, self.nw_file)
        report = """Input data:
  {}

Parameters:
  Minimum number of OTUs:\t{}
  Minimum sequence length:\t{}
  Long branch threshold:\t{}
  Minimum support value:\t{}
  Include:\t\t\t{}
  Exclude:\t\t\t{}
  Monophyly masking method:\t{}
  Rooting method:\t\t{}
  Outgroup rooting:\t\t{}
  Paralogy pruning method:\t{}
  Paralogy frequency threshold:\t{}
  Trim divergent percentage:\t{}
  Jackknife:\t\t\t{}""".format(input_files,
                               self.min_taxa,
                               self.min_len,
                               self.trim_lb,
                               self.min_support,
                               self.include,
                               self.exclude,
                               self.mask,
                               self.root,
                               self.outgroup,
                               self.prune,
                               self.trim_freq_paralogs,
                               self.trim_divergent,
                               self.jackknife)
        with open(log_path, "a") as log_file:
            log_file.write(report)
