"Module for storing and retrieving a set of options."

from phylopypruner import report

ON = u"[{}]".format(report.colorize("x", "green"))
OFF = u"[ ]"
BULLET = "~"
MISSING = "_"

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
        self._min_otu_occupancy = arguments.min_otu_occupancy
        self._min_gene_occupancy = arguments.min_gene_occupancy

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
        return self._include

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
    def min_otu_occupancy(self):
        "The mimum OTU occupancy allowed in the output."
        return self._min_otu_occupancy

    @min_otu_occupancy.setter
    def min_otu_occupancy(self, value):
        self._min_otu_occupancy = value

    @property
    def min_gene_occupancy(self):
        "The mimum amount of gene occupancy allowed in the output."
        return self._min_gene_occupancy

    @min_gene_occupancy.setter
    def min_gene_occupancy(self, value):
        self._min_gene_occupancy = value

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

    def to_str(self):
        "Return these settings as a formatted string."
        if self.min_len:
            len_str = self.min_len
        else:
            len_str = MISSING

        if self.min_support:
            support_str = int(self.min_support * 100)
        else:
            support_str = MISSING

        if self.min_otu_occupancy:
            min_otu_str = int(self.min_otu_occupancy * 100)
        else:
            min_otu_str = MISSING

        if self.min_gene_occupancy:
            min_gene_str = int(self.min_gene_occupancy * 100)
        else:
            min_gene_str = MISSING

        exclude_str = report.format_otus(self.exclude)
        outgroup_str = report.format_otus(self.outgroup)

        settings_report = u"""Parameters \
({} = required, {} = used, {} = unused):
  {} Always exclude the following OTUs: {}
  {} Remove sequences shorter than {} bases
  {} Remove branches longer than {} standard deviations of all branches
  {} Collapse branches with less support than {}% into polytomies
   {}  Mask monophylies using the {} method
  {} These OTUs were used for outgroup rooting: {}
  {} Alternative rooting method: {}
  {} The Paralogy frequency threshold is set to {}
  {} The Divergence threshold is set to {}
  {} Include these OTUs, even if deemed problematic: {}
   {}  Prune paralogs using the {} method
  {} Discard output alignments with fewer than {} sequences
  {} Taxon jackknifing is {}performed
  {} Discard OTUs with less than {}% occupancy
  {} Discard genes with less than {}% occupancy""".format(
      BULLET, ON, OFF,
      ON if self.exclude else OFF, exclude_str,
      ON if self.min_len else OFF, len_str,
      ON if self.trim_lb else OFF, self.trim_lb if self.trim_lb else MISSING,
      ON if self.min_support else OFF, support_str,
      BULLET, self.mask,
      ON if self.outgroup else OFF, outgroup_str,
      ON if self.root else OFF, self.root,
      ON if self.trim_freq_paralogs else OFF,
      self.trim_freq_paralogs if self.trim_freq_paralogs else MISSING,
      ON if self.trim_divergent else OFF,
      self.trim_divergent if self.trim_divergent else MISSING,
      ON if self.include else OFF, self.include,
      BULLET, self.prune,
      ON, self.min_taxa,
      ON if self.jackknife else OFF,
      "" if self.jackknife else "not ",
      ON if self.min_otu_occupancy else OFF, min_otu_str,
      ON if self.min_gene_occupancy else OFF, min_gene_str)
        return settings_report

    def print_settings(self):
        "Print these settings in a formatted report."
        print(self.to_str())
