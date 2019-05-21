"Module for storing and retrieving a set of options."

from phylopypruner import report

CHECKED = u"{}\u25CF{}".format(report.GREEN, report.NORMAL)
UNCHECKED = u"{}\u25CB{}".format(report.RED, report.NORMAL)
DEFAULT = "({}default{})".format(report.YELLOW, report.NORMAL)
# CHECKED = "{}{}[x]".format(report.GREEN, report.NORMAL)
# UNCHECKED = "{}{}[ ]".format(report.RED, report.NORMAL)
# DEFAULT = "({}default{})".format(report.YELLOW, report.NORMAL)

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
            len_str = "sequences shorter than {} bases were removed".format(self.min_len)
        else:
            len_str = "short sequences were not removed"
        if self.min_support:
            support_pct = int(self.min_support * 100)
            support_str = "collapse branches with less support than {}% into \
polytomies".format(support_pct)
        else:
            support_str = "do not collapse weakly supported branches into polytomies"

        settings_report = u"""
Parameters used/order of operations:
  \033[34m1. \033[0m {} {} {}
  \033[34m2. \033[0m {} minimum number of OTUs in output: {} {}
  \033[34m3. \033[0m {} branches longer than {} times the SD of all branches were removed
  \033[34m4. \033[0m {} {}
  \033[34m5. \033[0m {} Monophyly masking method.........: {}
  \033[34m6. \033[0m {} Rooting method...................: {}
  \033[34m7. \033[0m {} Outgroup rooting.................: {}
  \033[34m8. \033[0m {} the {} method was used for paralogy pruning {}
  \033[34m9. \033[0m {} The Paralogy frequency threshold.....: {}
  \033[34m10.\033[0m {} Trim divergent percentage........: {}
  \033[34m11.\033[0m {} include these OTUs, even if deemed problematic: {}
  \033[34m12.\033[0m {} always exclude the following OTUs: {}
  \033[34m13.\033[0m {} {}

Unused parameters:
  ...""".format(
      CHECKED if self.min_len else UNCHECKED, len_str, DEFAULT if not self.min_len
      else "",
      CHECKED if self.min_taxa else UNCHECKED, self.min_taxa, DEFAULT if
      self.min_taxa == 4 else "",
      CHECKED if self.trim_lb else UNCHECKED, self.trim_lb,
      CHECKED if self.min_support else UNCHECKED, support_str,
      CHECKED if self.mask else UNCHECKED, self.mask,
      CHECKED if self.root else UNCHECKED, self.root,
      CHECKED if self.outgroup else UNCHECKED, self.outgroup,
      CHECKED if self.prune else UNCHECKED, self.prune, DEFAULT if self.prune
      == "LS" else "",
      CHECKED if self.trim_freq_paralogs else UNCHECKED, self.trim_freq_paralogs,
      CHECKED if self.trim_divergent else UNCHECKED, self.trim_divergent,
      CHECKED if self.include else UNCHECKED, self.include,
      CHECKED if self.exclude else UNCHECKED, self.exclude,
      CHECKED if self.jackknife else UNCHECKED,
      "taxon jackknifing is performed" if self.jackknife else "do not perform \
taxon jackknifing")
        print(settings_report)
