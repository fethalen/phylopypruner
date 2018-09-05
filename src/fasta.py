"""
Module for working with the FASTA file format.
"""

from os import path
from multiple_sequence_alignment import MultipleSequenceAlignment

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
    msa = MultipleSequenceAlignment(filename, extension)

    with open(filename) as fasta_file:
        for line in fasta_file:
            line = line.rstrip()
            if line.startswith(">"):
                if sequence_data:
                    msa.add_sequence(None, description, sequence_data)
                description = line[1:]
                sequence_data = ""
            else:
                sequence_data += line
    return msa

def write(msa):
    """
    Takes a multiple sequence alignment object and a path as an input. Writes
    all sequences within that alignment to the provided filename.
    """
    with open(msa.filename, "w") as fasta_file:
        for sequence in msa.sequences:
            fasta_file.write(">{}\n".format(sequence.description))
            fasta_file.write("{}\n".format(sequence.sequence_data))
