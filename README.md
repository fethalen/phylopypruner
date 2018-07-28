PhyloPyPruner
-------------

PhyloPyPruner is a contamination-aware phylogenetic tree-based program for
orthology inference.

PhyloPyPruner is under active development. I appreciate if you try it on your
own data and [leave feedback](mailto:felix.thalen.1430@student.lu.se).

### Input data

PhyloPyPruner requires a multiple sequence alignment (MSA) and a Newick tree as
an input.

```
./phylopypruner [msa] [tree]
```

FASTA descriptions and Newick names should be in one of the following formats:
`OTU|ID` or `OTU@ID`, where `OTU` is the operational taxonomical unit (usually
a taxon name) and `ID` is a unique annotation or sequence identifier. For
example: `>Meiomenia_swedmarki|Contig00001_Hsp90`.

Sequence descriptions and tree names are not allowed to deviate from each
other. Sequence data needs to be [valid IUPAC nucleotide or amino acid
sequences](https://www.bioinformatics.org/sms/iupac.html).

© [Kocot Lab](https://www.kocotlab.com/) 2018
