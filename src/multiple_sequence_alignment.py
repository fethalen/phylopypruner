"""
Data type for storing a set of amino acid or nucleotide sequences.
"""

from sequence import Sequence

class MultipleSequenceAlignment(object):
    """
    Represents a set of sequences.
    """
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
        """
        The file extension used for this multiple sequence alignment.
        """
        return self._extension

    @extension.setter
    def extension(self, value):
        self._extension = value

    @property
    def sequences(self):
        """
        The sequences contained within this alignment.
        """
        return self._sequences

    @sequences.setter
    def sequences(self, value):
        self._sequences = value

    def add_sequence(self, sequence=None, description="", sequence_data=""):
        """
        Adds a new sequence to this alignment.
        """
        if not sequence:
            sequence = Sequence()

        sequence.description = description
        sequence.sequence_data = sequence_data

        self.sequences.append(sequence)
        return sequence

    def remove_sequence(self, sequence):
        """
        Remove the provided sequence from this alignment, if it exists.
        """
        try:
            self.sequences.remove(sequence)
        except ValueError:
            print("Sequence {} not found".format(sequence.description))

    def get_sequence(self, description):
        """
        Takes a FASTA description as an input and returns the matching sequence
        object, if a sequence with that description is found within this
        alignment.
        """
        for sequence in self.sequences:
            if description is sequence.description:
                return sequence
