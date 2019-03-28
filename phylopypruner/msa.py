"Data type for storing a set of amino acid or nucleotide sequences."

from __future__ import absolute_import
import re
from sys import stderr
from phylopypruner.sequence import Sequence

class MultipleSequenceAlignment(object):
    "Represents a set of sequences."
    def __init__(self, filename="", extension=None):
        self._filename = filename
        self._extension = extension
        self._sequences = []

    def __str__(self):
        return self.filename

    def __len__(self):
        return len(self.sequences)

    def __nonzero__(self):
        return True

    def __bool__(self):
        return True

    @property
    def filename(self):
        """The name of the file to which this multiple sequence alignment belongs.
        """
        return self._filename

    @filename.setter
    def filename(self, value):
        self._filename = value

    @property
    def extension(self):
        "The file extension used for this multiple sequence alignment."
        return self._extension

    @extension.setter
    def extension(self, value):
        self._extension = value

    @property
    def sequences(self):
        "Returns a list of all sequences in this alignment."
        return self._sequences

    @sequences.setter
    def sequences(self, value):
        self._sequences = value

    def add_sequence(self, sequence=None, description="", sequence_data=""):
        "Add a sequence object to the sequences list in this alignment object."
        if not sequence:
            sequence = Sequence()
            sequence.description = description
            sequence.sequence_data = sequence_data
            sequence.otu = re.split(r"\||@", sequence.description)[0]
            try:
                sequence.identifier = re.split(r"\||@", sequence.description)[1]
            except IndexError:
                print >> stderr, "# Warning: No description found on split with | or @"
                sequence.identifier = "None"
        elif description:
            sequence.otu = re.split(r"\||@", sequence.description)[0]
        elif sequence_data:
            sequence.identifier = re.split(r"\||@", sequence.description)[1]

        self.sequences.append(sequence)
        return sequence

    def get_sequence(self, description):
        """Takes a FASTA description as an input and returns the matching sequence
        object, if a sequence with that description is found within this
        alignment.
        """
        for sequence in self.sequences:
            if description == sequence.description:
                return sequence

    def iter_descriptions(self):
        """Returns an iterator object that includes all sequence descriptions
        in this alignment.
        """
        for sequence in self.sequences:
            yield sequence.description

    def iter_otus(self):
        "Returns an iterator object that includes all OTUs in this alignment."
        for sequence in self.sequences:
            yield sequence.otu

    def iter_identifiers(self):
        "Returns an iterator object that includes all IDs in this alignment."
        for sequence in self.sequences:
            yield sequence.identifier

    def gaps(self):
        """Returns the number of gap characters ('-', '?', or 'x') within this
        MultipleSequenceAlignment object.

        Returns
        -------
        gaps : int
            The number of gap characters within this MultipleSequenceAlignment
            object.
        """
        gaps = 0

        for sequence in self.sequences:
            gaps += self.alignment_len() - len(sequence.ungapped())

        return gaps

    def missing_data(self, otus_missing=0):
        """Returns the percent missing data within this
        MultipleSequenceAlignment object. The amount of missing data is
        calculated as follows: Sum the number of gap characters ('-', '?', or
        'x') within each alignment and then divide this number by the total
        number of positions within all alignments. Treat each OTU that is missing
        from the alignment (set by the user; DEFAULT: 0), as having gaps
        equivalent to the entire length of this MultipleSequenceAlignment
        object.

        Parameters
        ----------
        otus_missing : int
            The number of OTUs missing from this alignment (used to calculate
            missing data within a supermatrix; 0 by default).

        Returns
        -------
        missing_data : float
            The percent missing data within this alignment.
        """
        gaps = 0
        no_of_sequences = float(len(self))
        alignment_len = self.alignment_len()
        gaps += self.gaps()
        gaps += otus_missing * alignment_len

        if no_of_sequences > 0:
            return float(gaps) / (float(alignment_len) *
                                  float(no_of_sequences + otus_missing))
        else:
            return 0

    def otus(self):
        "Returns a list of the OTUs present within this alignment."
        otus_in_alignment = set()
        for sequence in self.sequences:
            otus_in_alignment.add(sequence.otu)
        return otus_in_alignment

    def alignment_len(self):
        """Returns the length of this alignment.

        Returns
        -------
        length : int
            This alignment's length.
        """
        return len(self.sequences[0])

    def is_dna(self):
        """Returns True if MSA consists of nucleotides (as opposed to amino
        acids). If either the character B or U is found, then the sequence must
        be DNA data, since these are not found in protein data. If these
        characters aren't encountered, however, then DNA/protein will be
        inferred by considering the amount of A, C, G, and T characters. If any
        sequence within this MSA contains at least 50% of these characters,
        then nucleotides are assumed. Please not that this method is not 100%
        accurate, especially for very short sequences.

        Returns
        -------
        nucleotides : bool
            True if this MSA consists of nucleotides, else False.
        """
        all_seqs = ""

        for sequence in self.sequences:
            all_seqs += sequence.sequence_data.lower()

        if "b" in all_seqs or "u" in all_seqs:
            return True

        all_seqs = all_seqs.replace("-", "")
        no_of_acgt = 0
        no_of_acgt += all_seqs.count("a")
        no_of_acgt += all_seqs.count("c")
        no_of_acgt += all_seqs.count("g")
        no_of_acgt += all_seqs.count("t")

        return no_of_acgt / len(all_seqs) >= 50
