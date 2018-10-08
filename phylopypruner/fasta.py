"""
Module for working with the FASTA file format.
"""

from __future__ import absolute_import
from os import path
from textwrap import wrap
from phylopypruner import msa

def read(filename):
    """
    Takes the path to a file in FASTA format and returns a list of Sequence
    objects for the sequences contained within that file. You can provide a
    data type (DNA/RNA/Protein) for the sequence or the data type will be
    guessed. Providing a data type can potentially speed up the process.
    """
    sequence_data = ""
    description = ""
    extension = path.splitext(filename)[1]
    alignment = msa.MultipleSequenceAlignment(filename, extension)

    with open(filename) as fasta_file:
        for line in fasta_file:
            line = line.rstrip()
            if line.startswith(">"):
                if sequence_data:
                    alignment.add_sequence(None, description, sequence_data)
                description = line[1:]
                sequence_data = ""
            else:
                sequence_data += line
        if sequence_data:
            alignment.add_sequence(None, description, sequence_data)
    return alignment

def write(alignment, max_column):
    """
    Takes a multiple sequence alignment object and a path as an input. Writes
    all sequences within that alignment to the provided filename.
    """
    with open(alignment.filename, "w") as fasta_file:
        for sequence in alignment.sequences:
            if max_column:
                seq_data = "\n".join(wrap(sequence.sequence_data, max_column))
            else:
                seq_data = sequence.sequence_data
            fasta_file.write(">{}\n".format(sequence.description))
            fasta_file.write("{}\n".format(seq_data))
