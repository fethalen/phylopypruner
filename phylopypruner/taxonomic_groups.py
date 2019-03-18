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
    """Takes the path to a file that contains one, or more, 'taxonomic groups'
    as an input. Each group is contained within a single line and consists of
    two parts: a group name and one, or more, taxa. Each line starts with the
    group name, followed by a colon (':'). Each taxa that follows is separated
    by a semi colon (';'). For example:

    Mollusca:Crepidula;Nucella;Aplysia;Phenacolepa
    Annelida:Helobdella;Capitella;Pectinaria;Clymenella
    Entoprocta:Loxosomella;Loxosoma;Barentsia;Pedicellina

    This format ensures that spaces (' ') are also supported in taxon names.

    Returns a list of TaxonomicGroup objects.
    """
    groups = []
    otus = []
    group = None

    with open(path) as groups_file:
        for line in groups_file:
            line = line.rstrip()

            try:
                entries = line.split(',')
                group_name = entries[0]
                otus = entries[1:]
            except:
                raise ValueError("couldn't parse taxonomic group file")

            group = TaxonomicGroup(group_name)
            group.otus = otus
            groups.append(group)

    return groups
