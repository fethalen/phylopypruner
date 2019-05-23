"Module for storing and retrieving a set of options."

from phylopypruner import report

ON = u"{}\u2713{}".format(report.GREEN, report.NORMAL)
OFF = u"{}x{}".format(report.RED, report.NORMAL)
DEFAULT = "({}default{})".format(report.YELLOW, report.NORMAL)
MISSING = report.underline("   ")

class Settings(object):
    """
    A list of settings used in a single run.
    """
    def __init__(self, arguments=None):
        self._fasta_file = None
        self._nw_file = None
        self._min_taxa = arguments.min_taxa
        self._min_len = arguments.min_len
        self._min_support = arguments.min_support
        self._trim_lb = arguments.trim_lb
        self._trim_zero_len = arguments.min_pdist
        self._outgroup = arguments.outgroup
        self._include = arguments.include
        self._exclude = arguments.exclude
        self._force_inclusion = arguments.force_inclusion
        self._root = arguments.root
        self._mask = arguments.mask
        self._prune = arguments.prune
        self._trim_freq_paralogs = arguments.trim_freq_paralogs
        self._trim_divergent = arguments.trim_divergent
        self._jackknife = arguments.jackknife
        self._taxonomic_groups = arguments.subclades

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

    @property
    def taxonomic_groups(self):
        "The path to a 'taxonomic groups' file."
        return self._taxonomic_groups

    @taxonomic_groups.setter
    def taxonomic_groups(self, value):
        self._taxonomic_groups = value

    @property
    def force_inclusion(self):
        """A list of OTUs, don't output orthologs where these OTUs are not
        present.
        """
        return self._force_inclusion

    @force_inclusion.setter
    def force_inclusion(self, value):
        self._force_inclusion = value

    def report(self, directory, log_path):
        "Generate a report of what settings were used."
        input_files = "Directory:\t{}".format(directory)
        report = """Input data:
  {}

Parameters used:
  Minimum number of OTUs...........: {}
  Minimum sequence length..........: {}
  Long branch threshold............: {}
  Minimum support value............: {}
  Include..........................: {}
  Exclude..........................: {}
  Monophyly masking method.........: {}
  Rooting method...................: {}
  Outgroup rooting.................: {}
  Paralogy pruning method..........: {}
  Paralogy frequency threshold.....: {}
  Trim divergent percentage........: {}
  Jackknife........................: {}""".format(
      input_files,
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
            return report

    def print_settings(self):
        """TODO: write in past or present?
        """
        if self.min_len:
            len_str = report.colorize(self.min_len, "cyan")
        else:
            len_str = MISSING
        if self.min_support:
            support_str = int(self.min_support * 100)
        else:
            support_str = MISSING

        settings_report = u"""Parameters \
({} = used, {} = unused, SDs = standard deviations):
{} remove sequences shorter than {} bases {}
{} remove branches longer than {} SDs of all branches {}
{} collapse branches with less support than {}% into polytomies {}
{} if 2+ sequences from a single OTU form a clade, keep the {} {}
{} OTUs used for outgroup rooting: {}
{} rooting used if outgroup rooting was off or failed: {}
{} the {} method was used for paralogy pruning {}
{} The Paralogy frequency threshold.....: {}
{} Trim divergent percentage........: {} {}
{} these OTUs were included, even if deemed problematic: {}
{} always exclude the following OTUs: {}
{} output alignments with fewer than {} sequences were discarded {}
{} taxon jackknifing is {}performed {}""".format(
    ON, OFF,
    ON if self.min_len else OFF, len_str, DEFAULT if not self.min_len
    else "",
    ON if self.trim_lb else OFF, self.trim_lb if self.trim_lb else MISSING,
    DEFAULT if not self.trim_lb else "",
    ON if self.min_support else OFF, support_str, DEFAULT if not
    self.min_support else "",
    ON,
    "sequence with the \n      shortest pairwise distance to its sister" if
    self.mask == "pdist" else "longest sequence",
    DEFAULT if self.mask == "pdist" else "",
    ON if self.outgroup else OFF, self.outgroup,
    ON if self.root else OFF, self.root,
    ON, self.prune, DEFAULT if self.prune
    == "LS" else "",
    ON if self.trim_freq_paralogs else OFF, self.trim_freq_paralogs,
    ON if self.trim_divergent else OFF, self.trim_divergent, DEFAULT if not
    self.trim_divergent else "",
    ON if self.include else OFF, self.include,
    ON if self.exclude else OFF, self.exclude,
    ON, self.min_taxa, DEFAULT if
    self.min_taxa == 4 else "",
    ON if self.jackknife else OFF,
    "" if self.jackknife else "not ",
    DEFAULT if not self.jackknife else "")
        print(settings_report)
