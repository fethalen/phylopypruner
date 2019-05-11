"Module for working with collections of Log objects."

from __future__ import absolute_import
import os, sys
import datetime
from phylopypruner import fasta
from textwrap import wrap
from collections import defaultdict
try:
    import matplotlib
    matplotlib.use("agg")
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import figure as fig
    import matplotlib.pylab as pylab
    MATPLOTLIB = True
except ImportError:
    print("{}tip:{} install Matplotlib (https://matplotlib.org/) to \
get a barplot of the paralog frequency".format("\033[92m", "\033[0m"))
    MATPLOTLIB = False
HAVE_DISPLAY = "DISPLAY" in os.environ
if not HAVE_DISPLAY:
    print("{}warning:{} no display found; can't generate plot".format(
        "\033[35m", "\033[0m"))
    MATPLOTLIB = False
TIMESTAMP = datetime.datetime.now().strftime("%Y-%m-%d")
SUM_HEADER = "id;alignments;sequences;otus;meanSequences;meanOtus;meanSeqLen;\
shortestSeq;longestSeq;pctMissingData;catAlignmentLen\n"
SUM_PATH = "/supermatrix_stats.csv"
OCCUPANCY_PLOT_FILE ="/occupancy_matrix.png"
FREQ_PLOT_FILE = "/paralogy_freq_plot.png"
FREQ_CSV_FILE = "/otu_stats.csv"
GAP_CHARACTERS = {"-", "?", "x"}

class Summary(object):
    "Represents a collection of Log objects from previous runs."
    def __init__(self, dir_out=None):
        self._logs = []

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
                        if not position in GAP_CHARACTERS:
                            presence[index] = True

                # Remove columns where no residues are present.
                for sequence in msa.sequences:
                    for index, position in enumerate(presence):
                        if not position:
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
            otus.update(log.msa.otus())

        no_of_otus = len(otus)

        for log in self.logs:
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
                    occupancy = len(sequence.ungapped()) / float(len(sequence))
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
        if MATPLOTLIB and not no_plot:
            plot_occupancy_matrix(occupancies, gene_partitions, otus, dir_out,
                                  below_threshold)

        return above_otu_threshold, above_gene_threshold

    def coverage_by_site(self):
        """...
        """
        for log in self.logs:
            for msa in log.msas_out:
                for sequence in msa.sequences:
                    print(len(sequence.ungapped()))

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
        if os.path.isfile(dir_out + FREQ_PLOT_FILE):
            os.remove(dir_out + FREQ_PLOT_FILE)

        if os.path.isfile(dir_out + FREQ_CSV_FILE):
            os.remove(dir_out + FREQ_CSV_FILE)

        threshold = None
        seen = set()
        paralog_freq = defaultdict(int) # key is OTU, value is no. of paralogs
        presence = defaultdict(int) # key is OTU, value is
        divergent = defaultdict(int) # key is OTU, divergent sequences is value
        first_iteration = True
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
            csv_out.write("otu;paralogyFrequency;timesAboveDivergenceThreshold\n")
            for otu, freq in paralog_freq.items():
                csv_out.write("{};{};{}\n".format(
                    otu, freq, divergent[otu]))

        otus = list(paralog_freq.keys())
        indexes = range(len(otus))
        freq = list(paralog_freq.values())

        if not MATPLOTLIB or no_plot:
            return paralog_freq

        try:
            plt.barh(y=indexes, color="black", edgecolor="black", width=freq, alpha=0.5)
        except TypeError:
            print >> sys.stderr, "# Plotting function lacking indices. Try --no-plot\n", indexes[:2]
            return paralog_freq
        plt.yticks(list(indexes), otus)
        plt.ylabel("OTU")
        plt.xlabel("number of paralogs / number of alignments OTU is in")
        plt.title("Paralogy Frequency")
        if threshold:
            plt.axvline(x=threshold,
                        color="red",
                        label="cutoff = {}".format(round(threshold, 3)),
                        linestyle="--")
            plt.legend(loc='upper right', fontsize=14)
        plot_figure = plt.gcf()
        plot_figure.set_size_inches(12.0, len(otus) * 0.17)
        plt.savefig(dir_out + FREQ_PLOT_FILE, bbox_inches='tight', dpi=300)

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
        homolog_report, row, stats = self.homolog_alignment_stats("Input Alignments", "input")

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
            A CSV-row where each column is separated by a semicolon, ';'.
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

        report = """
{}
{}
# of alignments..................:{:10d}
# of sequences...................:{:10d}
# of OTUs........................:{:10d}
avg # of sequences per alignment.:{:10d}
avg # of OTUs....................:{:10d}
avg sequence length (ungapped)...:{:10d}
shortest sequence (ungapped).....:{:10d}
longest sequence (ungapped)......:{:10d}
% missing data...................:{:10.2f}
concatenated alignment length....:{:10d}""".format(
    name,
    "-" * len(name),
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

        row = "{};{};{};{};{};{};{};{};{};{};{}\n".format(
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

        return report, row, stats

    def alignment_stats(self, name, title, homolog_stats):
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
            A CSV-row where each column is separated by a semicolon, ';'.
        """
        no_of_alignments = 0
        no_of_seqs = 0
        seq_lens = 0
        cat_alignment_len = 0
        pct_missing = 0.0
        shortest = None
        longest = None
        otus = set()

        for log in self.logs:
            for msa in log.msas_out:
                no_of_alignments += 1
                otus.update(msa.otus())
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
            for msa in log.msas_out:
                otus_missing = no_of_otus - len(list(msa.otus()))
                pct_missing += msa.missing_data(otus_missing)

        if no_of_alignments > 0:
            avg_no_of_seqs = int(no_of_seqs / no_of_alignments)
            missing_data = round((pct_missing / no_of_alignments) * 100, 1)
        else:
            avg_no_of_seqs = 0
            missing_data = 0

        if no_of_seqs > 0:
            avg_seq_len = int(seq_lens / no_of_seqs)
        else:
            avg_seq_len = 0

        report = """
{}
{}
# of alignments...................:{:10d} | {}
# of sequences....................:{:10d} | {}
# of OTUs.........................:{:10d} | {}
avg # of sequences per alignment..:{:10d} | {}
avg # of OTUs.....................:{:10d} | {}
avg sequence length (ungapped)....:{:10d} | {}
shortest sequence (ungapped)......:{:10d} | {}
longest sequence (ungapped).......:{:10d} | {}
% missing data....................:{:10.2f} | {}
concatenated alignment length.....:{:10d} | {}""".format(
    name,
    "-" * len(name) + " " * 29 + "  Input | Output",
    homolog_stats[1],
    no_of_alignments,
    homolog_stats[2],
    no_of_seqs,
    homolog_stats[3],
    no_of_otus,
    homolog_stats[4],
    avg_no_of_seqs,
    homolog_stats[5],
    avg_no_of_seqs,
    homolog_stats[6],
    avg_seq_len,
    homolog_stats[7],
    shortest,
    homolog_stats[8],
    longest,
    homolog_stats[9],
    missing_data,
    homolog_stats[10],
    cat_alignment_len)

        row = "{};{};{};{};{};{};{};{};{};{};{}\n".format(
            title,
            no_of_alignments,
            no_of_seqs,
            no_of_otus,
            avg_no_of_seqs,
            avg_no_of_seqs,
            avg_seq_len,
            shortest,
            longest,
            missing_data,
            cat_alignment_len)

        return report, row

    def report(self, title, dir_out, homolog_stats):
        """Output a summary of the files for this run.

        Parameters
        ----------
        title : str
            The ID of the summary file.
        dir_out : str
            Path to the output directory that the summary file is saved to.

        Return
        ------
        report : str
            Overview statistics of the summary.
        """
        report, row = self.alignment_stats("Alignments", title, homolog_stats)

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

def plot_occupancy_matrix(matrix, xlabels, ylabels, dir_out, below_threshold):
    """...
    """
    print("{}==>{} generating occupancy plot (disable with '--no-plot')".format(
        "\033[34m", "\033[0m"))

    plot_figure = plt.figure()
    axes = plot_figure.add_subplot(111)
    plot = axes.matshow(matrix, cmap="viridis_r", interpolation="nearest")
    plot_figure.colorbar(plot)

    # set default plot size and font size
    width = 8
    height = 8
    xfont = 7
    yfont = 7

    # allow for occupancy plots of various sizes
    if len(xlabels) > 100 or len(ylabels) > 100:
        width = 0.05 * len(xlabels)
        xfont = 3
        height = 0.05 * len(ylabels)
        yfont = 3

    axes.set_title("Occupancy Matrix")
    axes.xaxis.set_ticks_position("bottom")
    axes.set_xticks(list(range(len(xlabels))))
    axes.set_yticks(list(range(len(ylabels))))
    axes.set_xticklabels(list(xlabels), rotation="vertical",
                         fontsize=xfont, stretch="expanded")
    axes.set_yticklabels(list(ylabels), fontsize=yfont)

    otus_below, genes_below = below_threshold

    # Highlight gene partitions below the allowed threshold in red.
    for index in range(genes_below)[::-1]:
        axes.get_xticklabels()[-index - 1].set_color("red")

    # Highlight OTUs below the allowed threshold in red.
    for index in range(otus_below)[::-1]:
        axes.get_yticklabels()[-index - 1].set_color("red")

    plot_figure.set_size_inches(width, height)

    plot_figure.tight_layout()
    plt.xlabel("Gene partitions")
    plt.ylabel("OTUs")
    # Pad margins so that markers don't get clipped by the axes.
    plt.margins(0.2)
    # Tweak spacing to prevent clipping of tick-labels.
    plt.subplots_adjust(bottom=0.15)
    plt.savefig(dir_out + OCCUPANCY_PLOT_FILE, bbox_inches="tight", dpi=300)

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
