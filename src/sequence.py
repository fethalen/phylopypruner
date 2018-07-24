"""
Amino acid or nucleotide sequence data type.
"""

import re

IUPAC_CODES = {"nucleotide": re.compile("^[ATKMBVCNSWDGUYRH-.]*$"),
               "protein": re.compile("^[ABCDEFGHIKLMNPQRSTUVWYZX*-.]*$")}

class Sequence(object):
    """
    Represents a biological sequence. If no data type is provided, it will be
    determined based on the file's extension or the sequence content.
    """
    def __init__(self, description="", sequence_data="", data_type=None,
                 extension=None):
        self._description = str(description)
        self._sequence_data = str(sequence_data)
        self._data_type = self._guess_data_type(str(data_type))
        self._extension = str(extension)
        self._is_alignment = True if self.is_alignment else False

    def __str__(self):
        return self.description

    def __len__(self):
        return len(self.sequence_data)

    def __nonzero__(self):
        return True

    def __bool__(self):
        return True

    @property
    def description(self):
        """A description or name of the sequence."""
        return self._description

    @description.setter
    def description(self, value):
        self._description = value

    @property
    def sequence_data(self):
        """The raw sequence data."""
        return self._sequence_data

    @sequence_data.setter
    def sequence_data(self, value):
        self._sequence_data = value

    @property
    def data_type(self):
        """The type of the sequence (nucleotide/protein)."""
        return self._data_type

    @data_type.setter
    def data_type(self, value):
        if "DNA" or "RNA" or "Protein" in value:
            self._data_type = value
        else:
            raise ValueError("Data type needs to be nucleotide or protein.")

    @property
    def extension(self):
        """The sequence's file extension."""
        return self._extension

    @extension.setter
    def extension(self, value):
        self.extension = value

    @property
    def is_alignment(self):
        """
        Returns True if the sequence contain at least one gap character. Only
        dashes, '-' or dots '.' are considered to be gap characters.
        """
        return bool("-" in self.sequence_data)

    @is_alignment.setter
    def is_alignment(self, value):
        self._is_alignment = value

    def count(self, letter):
        """Returns the frequency of the provided letter in this sequence."""
        return self._sequence_data.count(letter)

    def _validate_sequence(self):
        if not self.sequence_data:
            return
        for letter in self.sequence_data:
            # if letter not in
            pass

    def _guess_data_type(self, data_type):
        """Determine the sequence's type."""
        if data_type:
            return data_type
        elif self.extension:
            if "fna" or "ffn" in self.extension.lower():
                if "U" in self.sequence_data:
                    return "RNA"
                else:
                    return "DNA"
            elif "faa" in self.extension.lower():
                data_type = "Protein"
        else:
            pass

    def gc_content(self):
        """
        Return this sequence's GC content as a floating point number.
        """
        sequence = self.sequence_data
        nucleotides = 'acgt'
        for base in seq_lower:
            if not base in nucleotides:
                print('Expected a nucleotide base (a, c, g, t), but found {}: \
                        '.format(base))
        gc_count = seq_lower.count('g') + seq_lower.count('c')
        return round(float(gc_count) / len(seq), 2)
