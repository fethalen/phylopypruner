"Module for working with collections of Log objects."

from __future__ import absolute_import
import os
import datetime
from phylopypruner import fasta
from textwrap import wrap
from collections import defaultdict
try:
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import figure as fig
    MATPLOTLIB = True
except ImportError:
    print("{}tip:{} install Matplotlib (https://matplotlib.org/) to \
get a barplot of the paralog frequency".format("\033[92m", "\033[0m"))
    MATPLOTLIB = False
HAVE_DISPLAY = "DISPLAY" in os.environ
if not HAVE_DISPLAY:
    print("{}warning:{} no display found; can't generate plot{}".format(
        "\033[35m", "\033[37m", "\033[0m"))
    MATPLOTLIB = False
TIMESTAMP = datetime.datetime.now().strftime("%Y-%m-%d")
SUM_HEADER = "id;alignments;sequences;otus;meanSequences;meanOtus;meanSeqLen;\
shortestSeq;longestSeq;pctMissingData;catAlignmentLen\n"
SUM_PATH = "/supermatrix_stats.csv"
FREQ_PLOT_FILE = "/paralogy_freq_plot.png"
FREQ_CSV_FILE = "/otu_stats.csv"

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

    def paralogy_frequency(self, dir_out, factor=None):
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

        if MATPLOTLIB:
            plt.barh(y=indexes, width=freq, color="black", edgecolor="black", alpha=0.5)
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
            fig = plt.gcf()
            fig.set_size_inches(12.0, len(otus) * 0.17)
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
        homolog_report, row = self.homolog_alignment_stats("Input Alignments", "input")

        with open(dir_out + SUM_PATH, "a") as sum_out_file:
            sum_out_file.write(row)

        return homolog_report

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
        pct_missing = 0.0
        shortest = None
        longest = None
        otus = set()

        for log in self.logs:
            msa = log.msa
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
            msa = log.msa
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
# of alignments:\t\t\t{}
# of sequences:\t\t\t\t{}
# of OTUs:\t\t\t\t{}
avg # of sequences per alignment:\t{}
avg # of OTUs:\t\t\t\t{}
avg sequence length (ungapped):\t\t{}
shortest sequence (ungapped):\t\t{}
longest sequence (ungapped):\t\t{}
% missing data:\t\t\t\t{}
concatenated alignment length:\t\t{}""".format(
    name,
    "-" * len(name),
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

    def alignment_stats(self, name, title):
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
# of alignments:\t\t\t{}
# of sequences:\t\t\t\t{}
# of OTUs:\t\t\t\t{}
avg # of sequences per alignment:\t{}
avg # of OTUs:\t\t\t\t{}
avg sequence length (ungapped):\t\t{}
shortest sequence (ungapped):\t\t{}
longest sequence (ungapped):\t\t{}
% missing data:\t\t\t\t{}
concatenated alignment length:\t\t{}""".format(
    name,
    "-" * len(name),
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

    def report(self, title, dir_out):
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
        report, row = self.alignment_stats("Output Alignments", title)

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
