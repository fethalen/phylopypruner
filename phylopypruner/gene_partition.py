"""
Module for working with gene partitions. Gene partitions are useful for
properly annotating Sequence object added to a Supermatrix object.
"""

class GenePartition(object):
    "Represents a gene partition."
    def __init__(self, name=None, start=None, stop=None):
        self._name = name if name else ""
        self._start = start if start else None
        self._stop = stop if stop else None

    def __str__(self):
        return "{} = {}-{}".format(self.name, self.start, self.stop)

    def __nonzero__(self):
        return True

    def __bool__(self):
        return True

    @property
    def name(self):
        "The name of this gene partition"
        return self._name

    @name.setter
    def name(self, value):
        self._name = value

    @property
    def start(self):
        "The start position of this gene partition."
        return self._start

    @start.setter
    def start(self, value):
        self._start = value

    @property
    def stop(self):
        "The stop position of this gene partition."
        return self._stop

    @stop.setter
    def stop(self, value):
        self._stop = value
