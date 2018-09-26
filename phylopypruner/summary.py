"A collection of Log objects from previous runs."

import os
import datetime
import fasta
from textwrap import wrap
from collections import defaultdict
try:
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import figure as fig
    MATPLOTLIB = True
except ImportError:
    print("{}suggestion{}: install Matplotlib (https://matplotlib.org/) to \
get a barplot of the paralog frequency".format("\033[92m\033[1m", "\033[0m"))
    MATPLOTLIB = False

TIMESTAMP = datetime.datetime.now().strftime("%Y-%m-%d")
SUM_HEADER = "id;alignments;sequences;sequences_avg;otus_avg;seq_len_avg;\
shortest_seq;longest_seq;pct_missing_data;cat_alignment_len\n"
SUM_PATH = "/{}_ppp_summary.csv".format(TIMESTAMP)
FREQ_PLOT_FILE = "/{}_ppp_paralog_freq.png".format(TIMESTAMP)
FREQ_CSV_FILE = "/{}_ppp_paralog_freq.csv".format(TIMESTAMP)

class Summary(object):
    """
    Represents a collection of Log objects from previous runs.
    """
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

    def paralogy_frequency(self, dir_out):
        if os.path.isfile(dir_out + FREQ_PLOT_FILE):
            os.remove(dir_out + FREQ_PLOT_FILE)

        if os.path.isfile(dir_out + FREQ_CSV_FILE):
            os.remove(dir_out + FREQ_CSV_FILE)

        seen = set()
        paralog_freq = defaultdict(int)
        for log in self.logs:
            for paralog in log.paralogs:
                otu = paralog.otu()
                # start counting at the first multiple of an OTU
                if otu in seen:
                    paralog_freq[otu] += 1
                seen.add(otu)

        with open(dir_out + FREQ_CSV_FILE, "w") as csv_out:
            csv_out.write("otu;paralogs\n")
            for otu, freq in paralog_freq.items():
                csv_out.write("{};{}\n".format(otu, freq))

        otus = paralog_freq.keys()
        freq = paralog_freq.values()

        if MATPLOTLIB:
            plt.barh(otus, freq, alpha=0.5)
            plt.ylabel("OTU")
            plt.xlabel("frequency")
            plt.title("Paralog Frequency")
            plt.tight_layout()
            fig = plt.gcf()
            fig.set_size_inches(18.5, 10.5)
            plt.savefig(dir_out + FREQ_PLOT_FILE, dpi=150)

        return paralog_freq

    def homolog_seq_files(self):
        "Returns the number of sequence files within the homolog phylogeny."
        no_of_files = 0
        no_of_files = sum(1 for log in self.logs)
        return no_of_files

    def homolog_seqs(self):
        "Returns the number of sequences within the homolog phylogeny."
        no_of_seqs = 0
        for log in self.logs:
            no_of_seqs += len(log.msa)
        return no_of_seqs

    def homolog_avg_seqs(self):
        "Returns the average number of sequences within the homolog phylogeny."
        seq_files = self.homolog_seq_files()
        if seq_files > 0:
            return self.homolog_seqs() / seq_files
        else:
            return 0

    def homolog_avg_seq_len(self):
        "Returns the average sequence length within the homolog phylogeny."
        seq_lens = 0
        sequences = 0
        for log in self.logs:
            for sequence in log.msa.sequences:
                sequences += 1
                seq_lens += len(sequence.ungapped())
        if sequences > 0:
            return seq_lens / sequences
        else:
            return 0

    def sequence_files(self):
        "Returns the number of sequence files contained within the Log objects."
        no_of_files = 0
        for log in self.logs:
            for ortholog in log.orthologs:
                if ortholog:
                    no_of_files += 1
        return no_of_files

    def sequences(self):
        "Returns the number of sequences within the Log objects."
        no_of_seqs = 0
        for log in self.logs:
            for msa_out in log.msas_out:
                no_of_seqs += len(msa_out)
        return no_of_seqs

    def shortest_sequence(self):
        "Returns the shortest sequence within all of the Log objects."
        shortest = None
        for log in self.logs:
            for msa_out in log.msas_out:
                for sequence in msa_out.sequences:
                    if not shortest or shortest > len(sequence.ungapped()):
                        shortest = len(sequence.ungapped())
        return shortest

    def longest_sequence(self):
        "Returns the longest sequence within all of the Log objects."
        longest = None
        for log in self.logs:
            for msa_out in log.msas_out:
                for sequence in msa_out.sequences:
                    if not longest or longest < len(sequence.ungapped()):
                        longest = len(sequence.ungapped())
        return longest

    def homolog_shortest_sequence(self):
        "Returns the shortest sequence within the homologs."
        shortest = None
        for log in self.logs:
            for sequence in log.msa.sequences:
                if not shortest or shortest > len(sequence.ungapped()):
                    shortest = len(sequence.ungapped())
        return shortest

    def homolog_longest_sequence(self):
        "Returns the longest sequence within all of the homologs."
        longest = None
        for log in self.logs:
            for sequence in log.msa.sequences:
                if not longest or longest < len(sequence):
                    longest = len(sequence.ungapped())
        return longest

    def avg_sequences(self):
        "Returns the average number of sequences per sequence file."
        no_of_files = self.sequence_files()
        if no_of_files > 0:
            return self.sequences() / self.sequence_files()
        else:
            return 0

    def avg_otus(self):
        "Returns the average number of OTUs in this summary."
        otus_total = 0
        ortholog_count = 0

        for log in self.logs:
            otus = None
            for ortholog in log.orthologs:
                if ortholog:
                    otus = set(ortholog.iter_otus())
                    otus_total += len(otus)
                    ortholog_count += 1

        if ortholog_count > 0:
            return otus_total / ortholog_count
        else:
            return 0

    def avg_seq_len(self):
        "Returns the average sequence length of all MSAs combined."
        seq_lens = 0
        sequences = 0
        for log in self.logs:
            for msa_out in log.msas_out:
                for sequence in msa_out.sequences:
                    sequences += 1
                    seq_lens += len(sequence.ungapped())
        if sequences > 0:
            return int(seq_lens / sequences)
        else:
            return 0

    def homolog_avg_otus(self):
        "Returns the average number of OTUs within the homologs."
        otus_total = 0
        ortholog_count = 0
        for log in self.logs:
            otus = set(log.tree.iter_otus())
            otus_total += len(otus)
            ortholog_count += 1
        if ortholog_count > 0:
            return otus_total / ortholog_count
        else:
            return 0

    def homolog_missing_data(self):
        "Returns the percent missing data within the homologs."
        sequences = 0
        pct_missing = 0.0
        for log in self.logs:
            for sequence in log.msa.sequences:
                sequences += 1
                ungapped = len(sequence.ungapped())
                missing = len(sequence) - ungapped
                if len(sequence) > 0:
                    pct_missing += float(missing) / float(len(sequence))
        if sequences > 0:
            return round(pct_missing / sequences, 3) * 100
        else:
            return 0

    def missing_data(self):
        "Returns the percent missing data within the orthologs."
        sequences = 0
        pct_missing = 0.0
        for log in self.logs:
            for msa_out in log.msas_out:
                for sequence in msa_out.sequences:
                    sequences += 1
                    ungapped = len(sequence.ungapped())
                    missing = len(sequence) - ungapped
                    if len(sequence) > 0:
                        pct_missing += float(missing) / float(len(sequence))
        if sequences > 0:
            return round(pct_missing / sequences, 3) * 100
        else:
            return 0

    def cat_alignment(self):
        """
        Returns the length of the concatenated alignment, if we were to
        concatenate all orthologs together.
        """
        seq_lens = 0
        for log in self.logs:
            for msa_out in log.msas_out:
                # It doesn't matter which sequence we picked since they're
                # aligned and missing positions are denoted by gaps.
                sequence = msa_out.sequences[0]
                seq_lens += len(sequence)
        return seq_lens

    def homolog_cat_alignment(self):
        """
        Returns the length of the concatenated alignment, if we were to
        concatenate all homologous sequences together.
        """
        seq_lens = 0
        for log in self.logs:
            sequence = log.msa.sequences[0]
            seq_lens += len(sequence)
        return seq_lens

    def homolog_report(self, dir_out):
        "Output summary statistics for the input files."
        homolog_report = """
Homolog report
--------------
# of alignments:\t\t\t{}
# of sequences:\t\t\t\t{}
avg # of sequences per alignment:\t{}
avg # of OTUs:\t\t\t\t{}
avg sequence length (ungapped):\t\t{}
shortest sequence (ungapped):\t\t{}
longest sequence (ungapped):\t\t{}
% missing data:\t\t\t\t{}
concatenated alignment length:\t\t{}""".format(
    self.homolog_seq_files(),
    self.homolog_seqs(),
    self.homolog_avg_seqs(),
    self.homolog_avg_otus(),
    self.homolog_avg_seq_len(),
    self.homolog_shortest_sequence(),
    self.homolog_longest_sequence(),
    self.homolog_missing_data(),
    self.homolog_cat_alignment())

        row = "{};{};{};{};{};{};{};{};{};{}\n".format(
            "homologs",
            self.homolog_seq_files(),
            self.homolog_seqs(),
            self.homolog_avg_seqs(),
            self.homolog_avg_otus(),
            self.homolog_avg_seq_len(),
            self.homolog_shortest_sequence(),
            self.homolog_longest_sequence(),
            self.homolog_missing_data(),
            self.homolog_cat_alignment())

        with open(dir_out + SUM_PATH, "a") as sum_out_file:
            sum_out_file.write(row)

        return homolog_report

    def report(self, title, dir_out):
        "Output a summary of the files for this run."
        report = """
Ortholog report
---------------
# of alignments:\t\t\t{}
# of sequences:\t\t\t\t{}
avg # of sequences per alignment:\t{}
avg # of OTUs:\t\t\t\t{}
avg sequence length (ungapped):\t\t{}
shortest sequence (ungapped):\t\t{}
longest sequence (ungapped):\t\t{}
% missing data:\t\t\t\t{}
concatenated alignment length:\t\t{}""".format(
    self.sequence_files(),
    self.sequences(),
    self.avg_sequences(),
    self.avg_otus(),
    self.avg_seq_len(),
    self.shortest_sequence(),
    self.longest_sequence(),
    self.missing_data(),
    self.cat_alignment())

        row = "{};{};{};{};{};{};{};{};{};{}\n".format(
            title,
            self.sequence_files(),
            self.sequences(),
            self.avg_sequences(),
            self.avg_otus(),
            self.avg_seq_len(),
            self.shortest_sequence(),
            self.longest_sequence(),
            self.missing_data(),
            self.cat_alignment())

        with open(dir_out + SUM_PATH, "a") as sum_out_file:
            sum_out_file.write(row)

        return report

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
