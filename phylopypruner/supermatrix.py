"Module for working with supermatrices."

from __future__ import absolute_import
import os
from phylopypruner.gene_partition import GenePartition
from phylopypruner.msa import MultipleSequenceAlignment
from phylopypruner.sequence import Sequence
from phylopypruner import fasta
SUPERMATRIX_FILENAME = "/supermatrix.fas"
PARTITION_FILENAME = "/partition_data.txt"

class Supermatrix(object):
    """
    Represents a supermatrix. A supermatrix consists of many multiple sequence
    alignments (MSAs), concatenated into a single alignment.
    """
    def __init__(self, dir_out):
        self._sequences = dict()
        self._gene_partitions = list()
        self._as_msa = MultipleSequenceAlignment(filename=dir_out +
                                                 SUPERMATRIX_FILENAME)

    def __len__(self):
        return len(self.sequences)

    def __nonzero__(self):
        return True

    def __bool__(self):
        return True

    @property
    def gene_partitions(self):
        """A set of GenePartition objects which belong to this Supermatrix object.
        """
        return self._gene_partitions

    @gene_partitions.setter
    def gene_partitions(self, value):
        self._gene_partitions = value

    @property
    def sequences(self):
        "A set of Sequence objects."
        return self._sequences

    @sequences.setter
    def sequences(self, value):
        self._sequences = value

    @property
    def as_msa(self):
        "Returns this supermatrix as a MultipleSequenceAlignment object."
        return self._as_msa

    @as_msa.setter
    def as_msa(self, value):
        self._as_msa = value

    def add_sequence(self, sequence):
        """Add a Sequence object as a gene partition to this Supermatrix
        object.

        Parameters
        ----------
        sequence : Sequence object
            The sequence you wish to add.

        Returns
        -------
        sequence : Sequence object
            The sequence that you added.
        """
        otu = sequence.otu

        if not otu in self.sequences:
            self.sequences[otu] = sequence.sequence_data
        else:
            self.sequences[otu] = self.sequences[otu] + sequence.sequence_data

        return sequence

    def add_partition(self, partition=None):
        """Add a GenePartition object to this Supermatrix object.

        Parameters
        ----------
        partition : GenePartition object
            The gene partition you wish to add.

        Returns
        -------
        partition : GenePartition object
            The gene partition that you added.
        """
        if not partition:
            partition = GenePartition()

        self.gene_partitions.append(partition)

    def partitions_from_summary(self, summary, dir_out, wrap=None):
        """Generate a set of partitions based on the output
        MultipleSequenceAlignment objects within the provided summary.

        Parameters
        ----------
        summary : Summary object
            The summary which contains the MultipleSequenceAlignment objects
            which you wish to add to this Supermatrix object.

        Returns
        -------
        self : Supermatrix
            This Supermatrix object with the MultipleSequenceAlignment objects
            added as gene partitions.
        """
        start = 1
        otus = set()

        # Infer whether the data consists of nucleotides or amino acids by
        # looking at the first log's input MSA.
        is_dna = summary.logs[0].msa.is_dna()

        for log in summary.logs:
            otus.update(log.msa.otus())

        for log in summary.logs:
            for msa in log.msas_out:
                gene_name = os.path.basename(msa.filename)
                gene_name = gene_name.split("_pruned")[0]
                stop = start + len(msa.sequences[0]) - 1
                self.add_partition(GenePartition(gene_name, start, stop))
                start = stop + 1

                otus_in_sequence = set()

                for sequence in msa.sequences:
                    self.add_sequence(sequence)
                    otus_in_sequence.add(sequence.otu)

                for otu in otus.difference(otus_in_sequence):
                    if is_dna:
                        filler = "N" * len(msa.sequences[0])
                    else:
                        filler = "X" * len(msa.sequences[0])
                    if not otu in self.sequences:
                        self.sequences[otu] = filler
                    else:
                        self.sequences[otu] = self.sequences[otu] + filler

        for otu in self.sequences:
            # Sequence expects the following format: OTU@identifier or
            # OTU|identifier, so we have to work around it.
            sequence = Sequence()
            sequence.description = otu
            sequence.otu = otu
            sequence.sequence_data = self.sequences[otu]
            self.as_msa.sequences.append(sequence)

        # Write the gene partitions to a file.
        with open(dir_out + PARTITION_FILENAME, "w") as partition_file:
            for gene_partition in self.gene_partitions:
                partition_file.write("AUTO, " + str(gene_partition) + "\n")

        fasta.write(self.as_msa, wrap)
