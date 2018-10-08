"""
A taxonomic group represents collections of operational taxonomical units
(OTUs).
"""

class TaxonomicGroup:
    """
    Represents a collection of operational taxonomical units (OTUs).
    """
    def __init__(self, name=""):
        self._name = name
        self._otus = []

    def __str__(self):
        return self.name

    def __len__(self):
        return len(self.otus)

    def __nonzero__(self):
        return True

    def __bool__(self):
        return True

    @property
    def name(self):
        "The name of this taxonomic group."
        return self._name

    @name.setter
    def name(self, value):
        self._name = value

    @property
    def otus(self):
        "A list of OTUs that belongs to this taxonomic group."
        return self._otus

    @otus.setter
    def otus(self, value):
        self._otus = value

    def add_otu(self, otu):
        "Add an OTU to this taxonomic group. Expects OTU to be a string."
        self.otus.append(otu)
        return otu

def read(path):
    """
    Takes a path to a file that describes one or more taxonomic groups as an
    input. Each group should begin with a single line that contains the group's
    name and end with a colon (':'). The groups name should be followed with a
    list of taxon names, separated by space (' '). Group and taxon names cannot
    contain any spaces.

    group_A:
    taxon_A1 taxon_A2 taxon_A3
    group_B:
    taxon_B1 taxon_B2 taxon_B3

    Returns a list of TaxonomicGroup objects.
    """
    groups = []
    otus = []
    group = None

    with open(path) as groups_file:
        for line in groups_file:
            line = line.rstrip()
            if line.endswith(":"):
                if otus:
                    group.otus = otus
                group_name = line.rstrip(":")
                group = TaxonomicGroup(group_name)
                groups.append(group)
                otus = []
            else:
                otus += [otu for otu in line.split(" ")]
        if otus:
            group.otus = otus

    return groups
