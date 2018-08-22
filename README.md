PhyloPyPruner
-------------

PhyloPyPruner is a contamination-aware phylogenetic tree-based program for
orthology inference.

PhyloPyPruner is under active development. I appreciate if you try it on your
own data and [leave feedback](mailto:felix.thalen.1430@student.lu.se).

### Input data

A multiple sequence alignment (MSA) and a Newick tree is required as an input.
The sequences in the MSA can be interleaved (span multiple lines).

```
./phylopypruner 16s.fas 16s.tre
```

FASTA descriptions and Newick names must match and has to be in one of the
following formats: `OTU|ID` or `OTU@ID`, where `OTU` is the operational
taxonomical unit (usually the species) and `ID` is a unique annotation or
sequence identifier. For example: `>Meiomenia_swedmarki|Contig00001_Hsp90`.

Sequence descriptions and tree names are not allowed to deviate from each
other. Sequence data needs to be [valid IUPAC nucleotide or amino acid
sequences](https://www.bioinformatics.org/sms/iupac.html).

© [Kocot Lab](https://www.kocotlab.com/) 2018
