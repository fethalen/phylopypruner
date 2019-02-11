"""
Module for working with supermatrices.
"""

class Supermatrix(object):
    """
    Represents a supermatrix. A supermatrix consists of many multiple sequence
    alignments (MSAs), concatenated into a single alignment.
    """
    def __init__(self):
        self._gene_partitions = list()

    def __len__(self):
        return len(self.sequence_data)

    def __nonzero__(self):
        return True

    def __bool__(self):
        return True

    @property
    def gene_partitions(self):
        """A set of GenePartition objects which belong to this Supermatrix object.
        """
        return self._gene_partitions

    @description.setter
    def description(self, value):
        self._description = value

    def add_sequence(self, sequence=None):
        """Add a Sequence object to this supermatrix.

        Parameters
        ----------
        sequence : Sequence
        """

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
            partition

        self.gene_partitions.append(partition)

    def from_genes
