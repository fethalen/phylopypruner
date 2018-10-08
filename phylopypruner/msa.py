"Data type for storing a set of amino acid or nucleotide sequences."

from __future__ import absolute_import
import re
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
        """
        The name of the file to which this multiple sequence alignment belongs.
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
            sequence.identifier = re.split(r"\||@", sequence.description)[1]
        elif description:
            sequence.otu = re.split(r"\||@", sequence.description)[0]
        elif sequence_data:
            sequence.identifier = re.split(r"\||@", sequence.description)[1]

        self.sequences.append(sequence)
        return sequence

    def get_sequence(self, description):
        """
        Takes a FASTA description as an input and returns the matching sequence
        object, if a sequence with that description is found within this
        alignment.
        """
        for sequence in self.sequences:
            if description == sequence.description:
                return sequence

    def iter_descriptions(self):
        """
        Returns an iterator object that includes all sequence descriptions in
        this alignment.
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

    def missing_data(self):
        """
        Returns the percentage of missing data of all sequences within this
        alignment.
        """
        no_of_sequences = len(self)
        missing_data = 0

        for sequence in self.sequences:
            missing_data += sequence.missing_data()

        if no_of_sequences > 0:
            return missing_data / no_of_sequences
        else:
            return 0

    def otus(self):
        "Returns a list of the OTUs present within this alignment."
        otus_in_alignment = set()
        for sequence in self.sequences:
            otus_in_alignment.add(sequence.otu)
        return otus_in_alignment

    def alignment_len(self):
        """
        Returns the length of this alignment.
        """
        return len(self.sequences[0])
