"Module for working with collections of Log objects."

from __future__ import absolute_import
import os
import datetime
from collections import defaultdict
from phylopypruner.report import underline
from phylopypruner import fasta
from phylopypruner import report
from phylopypruner import plot
try:
    import matplotlib as mpl
    if "DISPLAY" not in os.environ:
        mpl.use("agg")
    MATPLOTLIB = True
except ImportError:
    report.tip("install Matplotlib (https://matplotlib.org/) to generate \
plots")
    MATPLOTLIB = False
TIMESTAMP = datetime.datetime.now().strftime("%Y-%m-%d")
SUM_HEADER = "id,alignments,sequences,otus,meanSequences,meanOtus,meanSeqLen,\
shortestSeq,longestSeq,pctMissingData,catAlignmentLen\n"
SUM_PATH = "/supermatrix_stats.csv"
OCCUPANCY_PLOT_FILE = "/plots/occupancy_matrix.png"
FREQ_CSV_FILE = "/otu_stats.csv"
SPECIES_LIST_FILE = "/otu_list.txt"
GAP_CHARACTERS = {"-", "?", "x"}
RESOLUTION = 72
XFONT_SIZE = 2.5
YFONT_SIZE = 2.5

class Summary(object):
    "Represents a collection of Log objects from previous runs."
    def __init__(self):
        self._logs = []
        self._discarded_otus = []
        self._discarded_genes = []

    def __len__(self):
        return len(self.logs)

    def __nonzero__(self):
        return True

    def __bool__(self):
        return True

    @property
    def logs(self):
        "A list of Log objects."
        return self._logs

    @logs.setter
    def logs(self, value):
        self._logs = value

    @property
    def discarded_otus(self):
        """A set of OTUs that were discarded due to being below the allowed
        occupancy threshold.
        """
        return self._discarded_otus

    @discarded_otus.setter
    def discarded_otus(self, value):
        self._discarded_otus = value

    @property
    def discarded_genes(self):
        """A set of genes that were discarded due to being below the allowed
        gene occupancy threshold.
        """
        return self._discarded_genes

    @discarded_genes.setter
    def discarded_genes(self, value):
        self._discarded_genes = value

    def remove_gap_only_columns(self):
        """Iterate over the output alignments within this summary and remove
        positions where no residues are present.

        Parameters
        ----------
        None

        Returns
        -------
        self : Summary object
            This summary with the gap-only columns removed from the output
            MSAs.
        """
        for log in self.logs:
            for msa in log.msas_out:
                # Generate a list of booleans for each alignment: True = one
                # or more residues present at column, False = residues absent.
                presence = [False] * len(msa.sequences[0])

                # Identify columns where no residues are present.
                for sequence in msa.sequences:
                    for index, position in enumerate(sequence.sequence_data):
                        if position not in GAP_CHARACTERS:
                            presence[index] = True

                if False not in presence:
                    continue

                for sequence in msa.sequences:
                    for index in range(0, len(presence))[::-1]:
                        if not presence[index]:
                            sequence.sequence_data = remove_position_from_str(
                                sequence.sequence_data, index)

        return self

    def matrix_occupancy(self, dir_out, min_otu_occupancy=None,
                         min_gene_occupancy=None, no_plot=False):
        """Generate an occupancy matrix for each OTU and gene partition within
        this Summary object. Remove individual OTUs and gene partitions that
        are below the threshold, specified with 'min_otu_occupancy' and
        'min_gene_occupancy' respectively. Also call another function to plot
        the occupancy matrix and store it as a PNG file.

        Parameters
        ----------
        dir_out : str
            Store the PNG within this directory.
        min_otu_occupancy : float
            Minimum occupancy allowed per individual OTU. From 0.0 - 1.0
            (0-100%).
        min_gene_occupancy : float
            Minimum occupancy allowed per individual gene partition. From 0.0 -
            1.0 (0-100%).
        no_plot : bool
            Don't generate a plot of the occupancy matrix. Depending on the
            size of the supermatrix, this can save quite a lot of execution
            time.

        Returns
        -------
        above_otu_threshold : list
            A list of OTUs which have less occupancy than the user-provided
            threshold.
        above_gene_threshold : list
            A list of gene partitions which have less occupancy than the
            user-provided threshold.
        """
        occupancy_matrix = dict()
        otus = set()
        gene_partitions = list()
        gene_occupancies = dict()
        partitions_total = 0
        above_gene_threshold = list()

        for log in self.logs:
            for msa in log.msas_out:
                otus.update(msa.otus())

        no_of_otus = len(otus)

        for log in self.logs:
            if not log.msas_out:
                continue

            for msa in log.msas_out:

                gene_partition = os.path.basename(msa.filename)
                gene_partition = gene_partition.split("_pruned")[0]
                gene_partitions.append(gene_partition)
                otus_in_alignment = list()
                gene_occupancies[gene_partition] = 0
                partitions_total += 1

                # calculate the occupancy score for each present sequence
                for sequence in msa.sequences:
                    otu = sequence.otu
                    otus_in_alignment.append(otu)
                    if sequence:
                        occupancy = len(sequence.ungapped()) / float(len(sequence))
                    else:
                        occupancy = 0
                    gene_occupancies[gene_partition] += occupancy

                    if not otu in occupancy_matrix:
                        occupancy_matrix[otu] = [occupancy]
                    else:
                        occupancy_matrix[otu].append(occupancy)

                if min_gene_occupancy and no_of_otus > 0:
                    if (gene_occupancies[gene_partition] / no_of_otus) < min_gene_occupancy:
                        above_gene_threshold.append(msa)

                # add '0.0' for each missing sequence
                for otu in otus:
                    if otu in otus_in_alignment:
                        continue
                    if not otu in occupancy_matrix:
                        occupancy_matrix[otu] = [0.0]
                    else:
                        occupancy_matrix[otu].append(0.0)

        # generate a list of OTUs which have less occupancy than the threshold
        above_otu_threshold = list()
        if min_otu_occupancy and partitions_total > 0:
            for otu in occupancy_matrix:
                if sum(occupancy_matrix[otu]) / partitions_total < min_otu_occupancy:
                    above_otu_threshold.append(otu)

        # sort gene partitions by occupancy (from highest to lowest)
        sequence = list()
        for gene_partition, _ in sorted(gene_occupancies.items(), key=lambda item:
                                        (item[1], item[0]), reverse=True):
            sequence.append(gene_partitions.index(gene_partition))

        gene_partitions = [gene_partitions[index] for index in sequence]

        # sort OTUs by occupancy (from highest to lowest)
        occupancies = list()
        otus = list()

        for otu, occupancy in sorted(occupancy_matrix.items(), key=lambda item:
                                     (item[1], item[0]), reverse=True):
            otus.append(otu)
            occupancies.append(occupancy)

        # replace the old order of gene partitions with the new one
        for index, row in enumerate(occupancies):
            occupancies[index] = [row[column] for column in sequence]

        below_threshold = (len(above_otu_threshold), len(above_gene_threshold))

        # plot the occupancy matrix
        if MATPLOTLIB and not no_plot and occupancies:
            plot.occupancy_matrix(occupancies, gene_partitions, otus, dir_out,
                                  below_threshold)

        return above_otu_threshold, above_gene_threshold

    def paralogy_frequency(self, dir_out, factor=None, no_plot=False):
        """Calculate the paralogy frequency (PF) for all OTUs within this
        summary, where paralogy frequency is the number of paralogs divided by
        the number of alignments that each OTU is present in. Output the
        statistics in form of a CSV file, a PNG plot and return a dictionary.

        Parameters
        ----------
        dir_out : str
            Path to the directory to which you wish to output the PF plot.

        factor : float
            Used to visualize where the threshold for frequent paralogy cutoff
            was set.

        Returns
        -------
        paralog_freq : dict
            Dictionary where a OTU string is key and PF, as a float, is the value.
        """
        if os.path.isfile(dir_out + FREQ_CSV_FILE):
            os.remove(dir_out + FREQ_CSV_FILE)

        threshold = None
        seen = set()
        paralog_freq = defaultdict(int) # key is OTU, value is no. of paralogs
        presence = defaultdict(int) # key is OTU, value is
        divergent = defaultdict(int) # key is OTU, divergent sequences is value

        for log in self.logs:
            otus_in_alignment = log.msa.otus()
            for otu in otus_in_alignment:
                if otu:
                    presence[otu] += 1

            for paralog in log.paralogs:
                otu = paralog.otu()
                # start counting at the first multiple of an OTU
                if otu and otu in seen:
                    paralog_freq[otu] += 1
                seen.add(otu)

            for otu in log.divergent:
                divergent[otu] += 1

        # normalize paralogy frequency by how often the OTU is present
        for otu in presence:
            if not otu in paralog_freq:
                paralog_freq[otu] = 0
            else:
                paralog_freq[otu] = round(
                    (float(paralog_freq[otu]) / float(presence[otu])) * 100, 1)

        if factor:
            threshold = _std(list(paralog_freq.values())) * factor

        with open(dir_out + FREQ_CSV_FILE, "w") as csv_out:
            csv_out.write("otu,paralogyFrequency,timesAboveDivergenceThreshold\n")
            for otu, freq in paralog_freq.items():
                csv_out.write("{},{},{}\n".format(
                    otu, freq, divergent[otu]))

        # In addition to the OTU-stats, create a list of all the OTUs, which
        # can be useful for automating subsequent steps.
        with open(dir_out + SPECIES_LIST_FILE, "w") as species_list:
            for otu, _ in paralog_freq.items():
                species_list.write("{}\n".format(otu))

        otus = list(paralog_freq.keys())
        indices = range(len(otus))
        freq = list(paralog_freq.values())

        if MATPLOTLIB and not no_plot:
            plot.paralogy_frequency(indices, freq, otus, threshold, dir_out)

        return paralog_freq

    def homolog_report(self, dir_out):
        """Output a summary of the input alignments. Title is always
        'homologs'.

        Parameters
        ----------
        dir_out : str
            Path to the output directory that the summary file is saved to.

        Return
        ------
        report : str
            Overview statistics of the summary.
        """
        row, stats = self.homolog_alignment_stats("Input Alignments", "input")

        with open(dir_out + SUM_PATH, "a") as sum_out_file:
            sum_out_file.write(row)

        return stats

    def homolog_alignment_stats(self, name, title):
        """Returns a report, in text, and a CSV-row of a set of statistics
        based on this summary.

        Parameters
        ----------
        name : str
            The title to display in the text report.
        title : str
            The identifier for the CSV-row.

        Returns
        -------
        report : str
            A report of the statistics within this file.
        row : str
            A CSV-row where each column is separated by a semicolon, ','.
        """
        no_of_alignments = 0
        no_of_seqs = 0
        seq_lens = 0
        cat_alignment_len = 0
        otus_total = 0
        pct_missing = 0.0
        shortest = None
        longest = None
        otus = set()

        for log in self.logs:
            msa = log.msa
            no_of_alignments += 1
            otus.update(msa.otus())
            otus_total += len(msa.otus())

            cat_alignment_len += msa.alignment_len()

            for sequence in msa.sequences:
                no_of_seqs += 1
                seq_len = len(sequence.ungapped())
                seq_lens += seq_len

                if not shortest or shortest > seq_len:
                    shortest = seq_len

                if not longest or longest < seq_len:
                    longest = seq_len

        no_of_otus = len(otus)
        for log in self.logs:
            msa = log.msa
            otus_missing = no_of_otus - len(list(msa.otus()))
            pct_missing += msa.missing_data(otus_missing)

        if no_of_alignments > 0:
            avg_no_of_seqs = int(no_of_seqs / no_of_alignments)
            avg_no_of_otus = int(otus_total / no_of_alignments)
            missing_data = round((pct_missing / no_of_alignments) * 100, 1)
        else:
            avg_no_of_seqs = 0
            avg_no_of_otus = 0
            missing_data = 0

        if no_of_seqs > 0:
            avg_seq_len = int(seq_lens / no_of_seqs)
        else:
            avg_seq_len = 0

        if not shortest:
            shortest = 0
        if not longest:
            longest = 0

        row = "{},{},{},{},{},{},{},{},{},{},{}\n".format(
            title,
            no_of_alignments,
            no_of_seqs,
            no_of_otus,
            avg_no_of_seqs,
            avg_no_of_otus,
            avg_seq_len,
            shortest,
            longest,
            missing_data,
            cat_alignment_len)

        stats = [title, no_of_alignments, no_of_seqs, no_of_otus,
                 avg_no_of_seqs, avg_no_of_otus, avg_seq_len, shortest, longest,
                 missing_data, cat_alignment_len]

        return row, stats

    def alignment_stats(self, name, title, homolog_stats=None):
        """Returns a report, in text, and a CSV-row of a set of statistics
        based on this summary.

        Parameters
        ----------
        name : str
            The title to display in the text report.
        title : str
            The identifier for the CSV-row.
        homolog_stats : list
            Optional statistics for homlogs. If these are not provided, then no
            report is returned.

        Returns
        -------
        report : str
            A report of the statistics within this file.
        row : str
            A CSV-row where each column is separated by a semicolon, ','.
        """

        pct_missing = 0.0
        otus = set()
        stats = {
            "alignments": 0,
            "sequences": 0,
            "seq_lens": 0,
            "cat_alignment_len": 0,
            "short": 0,
            "long": 0,
            "ultrashort": 0,
            "divergent": 0,
            "shortest": 0,
            "longest": 0,
            "collapsed": 0,
            "avg_no_of_seqs": 0,
            "avg_seq_len": 0,
            "no_of_otus": 0,
            "missing_data": 0,
            "otus_removed_count": len(self.discarded_otus),
            "otus_removed_pct": 0,
            "genes_removed_count": len(self.discarded_genes),
            "genes_removed_pct": 0
                }

        for log in self.logs:
            for msa in log.msas_out:
                stats["alignments"] += 1
                otus.update(msa.otus())
                stats["cat_alignment_len"] += msa.alignment_len()

                for sequence in msa.sequences:
                    stats["sequences"] += 1
                    seq_len = len(sequence.ungapped())
                    stats["seq_lens"] += seq_len

                    if stats["shortest"] == 0 or stats["shortest"] > seq_len:
                        stats["shortest"] = seq_len

                    if stats["longest"] < seq_len:
                        stats["longest"] = seq_len

        stats["no_of_otus"] = len(otus)

        for log in self.logs:
            stats["short"] += len(log.trimmed_seqs)
            stats["long"] += len(log.lbs_removed)
            stats["ultrashort"] += len(log.ultra_short_branches)
            stats["divergent"] += len(log.divergent_removed)
            stats["collapsed"] += log.collapsed_nodes

            for msa in log.msas_out:
                otus_missing = stats["no_of_otus"] - len(list(msa.otus()))
                pct_missing += msa.missing_data(otus_missing)

        if stats["alignments"] > 0:
            stats["avg_no_of_seqs"] = int(
                stats["sequences"] / stats["alignments"])
            stats["missing_data"] = round(
                (pct_missing / stats["alignments"]) * 100, 1)

        if stats["sequences"] > 0:
            stats["avg_seq_len"] = int(stats["seq_lens"] / stats["sequences"])

        # determine the column width from the longest statistic
        if homolog_stats:
            col_width = max(
                len(str(stats["sequences"])),
                len(str(stats["alignments"])),
                len(str(stats["cat_alignment_len"])),
                len(str(stats["no_of_otus"])),
                len(str(stats["longest"])),
                len(str(homolog_stats[1])),
                len(str(homolog_stats[2])),
                len(str(homolog_stats[8])),
                len(str(homolog_stats[10]))
                ) + 1
        else:
            col_width = max(
                len(str(stats["sequences"])),
                len(str(stats["alignments"])),
                len(str(stats["cat_alignment_len"])),
                len(str(stats["no_of_otus"])),
                len(str(stats["longest"])),
                ) + 1

        if col_width < 4:
            col_width = 4

        if otus:
            stats["otus_removed_pct"] = 100 * stats["otus_removed_count"] / len(otus)

        header = "Alignment statistics:\n  " +\
                underline("{:33s}  {:{}s}   {:{}s}".format(
                    name, "Input", col_width, "Output", col_width))
        methods_header = "Methods summary:\n  " +\
                underline("{:29s} {:{}s}    {:{}s}".format(
                    name, "No. removed", col_width, "% of input", col_width))

        row = "{},{},{},{},{},{},{},{},{},{},{}\n".format(
            title,
            stats["alignments"],
            stats["sequences"],
            stats["no_of_otus"],
            stats["avg_no_of_seqs"],
            stats["avg_no_of_seqs"],
            stats["avg_seq_len"],
            stats["shortest"],
            stats["longest"],
            stats["missing_data"],
            stats["cat_alignment_len"])


        if not homolog_stats:
            return "", row

        stats_report = """
{}
  No. of alignments                 {:{}d}   {:{}d}
  No. of sequences                  {:{}d}   {:{}d}
  No. of OTUs                       {:{}d}   {:{}d}
  Avg no. of sequences / alignment  {:{}d}   {:{}d}
  Avg no. of OTUs / alignment       {:{}d}   {:{}d}
  Avg sequence length (ungapped)    {:{}d}   {:{}d}
  Shortest sequence (ungapped)      {:{}d}   {:{}d}
  Longest sequence (ungapped)       {:{}d}   {:{}d}
  % missing data                    {:{}.2f}   {:{}.2f}
  Concatenated alignment length     {:{}d}   {:{}d}

{}
  Short sequences            {:{}d}         {:{}.2f}
  Long branches              {:{}d}         {:{}.2f}
  Ultrashort distance pairs  {:{}d}         {:{}.2f}
  Divergent sequences        {:{}d}         {:{}.2f}
  Collapsed nodes            {:{}d}         {:{}.2f}
  OTUs < occupancy threshold {:{}d}         {:{}.2f}
  Genes < occupancy threshold{:{}d}         {:{}.2f}""".format(
      header,
      homolog_stats[1], col_width, stats["alignments"], col_width,
      homolog_stats[2], col_width, stats["sequences"], col_width,
      homolog_stats[3], col_width, stats["no_of_otus"], col_width,
      homolog_stats[4], col_width, stats["avg_no_of_seqs"], col_width,
      homolog_stats[5], col_width, stats["avg_no_of_seqs"], col_width,
      homolog_stats[6], col_width, stats["avg_seq_len"], col_width,
      homolog_stats[7], col_width, stats["shortest"], col_width,
      homolog_stats[8], col_width, stats["longest"], col_width,
      homolog_stats[9], col_width, stats["missing_data"], col_width,
      homolog_stats[10], col_width, stats["cat_alignment_len"], col_width,
      methods_header,
      stats["short"], col_width, 100 * stats["short"] / homolog_stats[2], col_width,
      stats["long"], col_width, 100 * stats["long"] / homolog_stats[2], col_width,
      stats["ultrashort"], col_width, 100 * stats["ultrashort"] / homolog_stats[2], col_width,
      stats["divergent"], col_width, 100 * stats["divergent"] / homolog_stats[2], col_width,
      stats["collapsed"], col_width, 100 * stats["collapsed"] / homolog_stats[2],
      col_width,
      stats["otus_removed_count"], col_width, stats["otus_removed_pct"], col_width,
      stats["genes_removed_count"], col_width,
      100 * stats["genes_removed_count"] / homolog_stats[2], col_width
  )

        return stats_report, row

    def report(self, title, dir_out, homolog_stats=None):
        """Output a summary of the files for this run.

        Parameters
        ----------
        title : str
            The ID of the summary file.
        dir_out : str
            Path to the output directory that the summary file is saved to.
        homolog_stats : str
            Provide this string if you wish to get a (printable) report.

        Return
        ------
        report : str
            Overview statistics of the summary.
        """
        if homolog_stats:
            report, row = self.alignment_stats("Description", title, homolog_stats)
        else:
            report, row = self.alignment_stats("Description", title)

        with open(dir_out + SUM_PATH, "a") as sum_out_file:
            sum_out_file.write(row)

        return report

    def otus(self):
        "Returns a set of all OTUs within this Summary object."
        otus_in_summary = set()
        for log in self.logs:
            for msa_out in log.msas_out:
                otus_in_summary.update(msa_out.otus())

        return otus_in_summary

    def write_msas(self, wrap=None):
        """
        Write the alignments within the Log objects in this summary's logs, if any,
        into one or more alignments file. If wrap has been specified, wrap
        sequence data at the provided column.
        """
        for log in self.logs:
            for msa in log.msas_out:
                fasta.write(msa, wrap)

def mk_sum_out_title(dir_out):
    """
    Takes the path to an output directory as an input and writes
    """
    if os.path.isfile(dir_out + SUM_PATH):
        os.remove(dir_out + SUM_PATH)

    with open(dir_out + SUM_PATH, "w") as sum_out_file:
        sum_out_file.write(SUM_HEADER)

def remove_position_from_str(string, position):
    """Takes a string and a position as an input and removes the characters at
    the provided position from the string.

    Parameters
    ----------
    string : str
        The string which you want to perform this action on.
    position : int
        The single position you which to remove.

    Returns
    -------
    string : str
        The input string with the provided position removed.
    """
    return string[:position] + string[position + 1:]

def _mean(data):
    """Returns the sample arithmetic mean of data. 0 is returned if an empty
    list was provided.

    Parameters
    ----------
    data : list of floats

    Returns
    _______
    out: float
        The sample arithmetic mean of data.
    """
    return float(sum(data)) / max(len(data), 1)

def _sdm(data):
    """Returns the squared deviations from the mean (SDM) of data.

    Parameters
    ----------
    data : list of floats

    Returns
    -------
    out : float
        The sum of square deviations of data.
    """
    return sum((x - _mean(data))**2 for x in data)

def _std(data):
    """Return the population standard deviation of data.

    Parameters
    ----------
    data : list of floats

    Returns
    -------
    out : float
        The population standard deviation of data.
    """
    if len(data) < 2:
        raise ValueError('variance requires at least two data points')
    return (_sdm(data) / len(data)) ** 0.5
