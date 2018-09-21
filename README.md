PhyloPyPruner
-------------

PhyloPyPruner is a contamination-aware phylogenetic tree-based program for
orthology inference.

PhyloPyPruner is under active development. I appreciate if you try it on your
own data and [leave feedback](mailto:felix.thalen.1430@student.lu.se).

To get a list of options, either run the software without any parameters or by
using the `-h` or `--help` flag.

### Input data

Either provide a single multiple sequence alignment (MSA) and a Newick tree by
using the `--msa` and `--tree` flags:

```
./phylopypruner --msa 16s.fas --tree 16s.tre
```

or, provide a `path` to an input directory, containing multiple trees and
alignments, by typing `--dir path`.

FASTA descriptions and Newick names must match and has to be in one of the
following formats: `OTU|ID` or `OTU@ID`, where `OTU` is the operational
taxonomical unit (usually the species) and `ID` is a unique annotation or
sequence identifier. For example: `>Meiomenia_swedmarki|Contig00001_Hsp90`.

Sequence descriptions and tree names are not allowed to deviate from each
other. Sequence data needs to be [valid IUPAC nucleotide or amino acid
sequences](https://www.bioinformatics.org/sms/iupac.html).

For inputting multiple files, you provide a path to the directory in which
these files reside.

```
./phylopypruner --dir <path>
```

The program will automatically look for trees and alignments with the same name
and run for each of these pair.

### Options

| **Command** | **Meaning** |
|----|----|
| `-h, --help` | display a help message and exit |
| `-v, --version` | display the version number and exit |
| `--msa` | path to a multiple sequence alignment (MSA) in FASTA format |
| `--tree` | path to a Newick tree file |
| `--dir` | path to a directory that contains multiple trees and alignments |
| `--output` | path to the output directory |
| `--min-taxa <threshold>` | minimum number of OTUs allowed in output |
| `--min-seq <threshold>` | remove sequences shorter than `<threshold>` |
| `--min-support <threshold>` | collapse nodes with a support value below `<threshold>` into polytomies |
| `--trim-lb <factor>` | remove sequences with a branch length that is `<factor>` times longer than the standard deviation of all branches |
| `--outgroup <OTU> [<OTU> ...]` | one or more outgroup OTUs; root using outgroup rooting, if the outgroup forms a non-repetetive, monophyletic group |
| `--root {midpoint,clock,MAD}` | reroot tree using this rooting method if an outgroup hasn't been provided or if outgroup rooting fails |
| `--mask {longest,pdist}` | mask monophylies using this method |
| `--prune {LS,MI,MO,RT,1to1}` | prune paralogs using this method |
| `--wrap <max column>` | wrap output sequences at column `<max column>` |

### Output files

The following files are generated after running this program:

* `<orthologs>/...` – output alignments
* `<timestamp>_ppp_run.log` – detailed report of each actions performed
* `<timestamp>_ppp_ortho_stats.csv` – statistics for output alignments
* `<timestamp>_ppp_paralog_freq.csv` – paralogy frequency data
* `<timestamp>_ppp_paralog_freq.png` – paralogy frequency plot\*

\* – only produced if [Matplotlib](https://matplotlib.org/) is installed

For a more details, please see [the
Wiki](https://gitlab.com/fethalen/phylopypruner/wikis).

© [Kocot Lab](https://www.kocotlab.com/) 2018
