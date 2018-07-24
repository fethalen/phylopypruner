"""
Module for working with the FASTA file format.
"""

from os import path
from sequence import Sequence

def read(filename, data_type=None):
    """
    Takes the path to a file in FASTA format and returns a list of Sequence
    objects for the sequences contained within that file. You can provide a
    data type (DNA/RNA/Protein) for the sequence or the data type will be
    guessed. Providing a data type can potentially speed up the process.
    """
    sequences = []
    sequence = ""
    description = ""
    data_type = data_type
    extension = path.splitext(filename)[1]

    with open(filename) as fasta_file:
        for line in fasta_file:
            line = line.rstrip()
            if line.startswith(">"):
                if sequence:
                    sequences.append(Sequence(description, sequence, data_type,
                                              extension))
                description = line
            else:
                sequence += line
    return sequences

def write(sequences, path):
    """
    Takes a list of sequences and a path as an input and writes that sequence
    data to a file in the provided path.
    """
    pass
