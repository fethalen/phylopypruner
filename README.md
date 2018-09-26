PhyloPyPruner
-------------

PhyloPyPruner is a phylogenetic tree-based orthology inference program for
refining orthology inference made by graph-based (or phenetic) approaches.

PhyloPyPruner is under active development. I appreciate if you try it on your
own data and [leave feedback](mailto:felix.thalen.1430@student.lu.se).

To get a list of options, either run the software without any parameters or by
using the `-h` or `--help` flag. For more details, please see [the
Wiki](https://gitlab.com/fethalen/phylopypruner/wikis).

### Features

* Remove short sequences
* Remove sequences with long branches
* Collapse weakly supported nodes into polytomies
* Five different paralogy pruning algorithms
* Measure and remove OTUs with frequent paralogs
* Get rid of problematic OTUs using taxon jackknifing
* Mask monophylies by choosing the longest sequence or using pairwise distance

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

### Output files

The following files are generated after running this program:

* `<timestamp>_<orthologs>/...` – output alignments
* `<timestamp>_ppp_summary.log` – summary statistics for all alignments
* `<timestamp>_ppp_run.log` – detailed report of each performed action
* `<timestamp>_ppp_ortho_stats.csv` – statistics for output alignments
* `<timestamp>_ppp_paralog_freq.csv` – paralogy frequency data
* `<timestamp>_ppp_paralog_freq.png` – paralogy frequency plot\*

\* – only produced if [Matplotlib](https://matplotlib.org/) is installed

© [Kocot Lab](https://www.kocotlab.com/) 2018
