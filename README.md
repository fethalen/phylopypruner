PhyloPyPruner
-------------

PhyloPyPruner is a tree-based orthology inference program for refining
orthology inference made by a graph-based approach. In addition to implementing
previously published paralogy pruning algorithms seen in
[PhyloTreePruner](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3825643/),
[UPhO](https://academic.oup.com/mbe/article/33/8/2117/2578877),
[Agalma](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3840672/) and
[Phylogenomic Dataset
Reconstruction](https://www.ncbi.nlm.nih.gov/pubmed/25158799), this software
also provides methods for identifying and getting rid of operational
taxonomical units (OTUs) that display contamination-like issues.

PhyloPyPruner is currently under active development and I would appreciate it
if you try this software on your own data and [leave
feedback](mailto:felix.thalen.1430@student.lu.se).

See [the Wiki](https://gitlab.com/fethalen/phylopypruner/wikis) for more
details.

## Features

* Remove short sequences
* Remove relatively long branches
* Collapse weakly supported nodes into polytomies
* Prune paralogs using one out of five methods
* Measure paralogy frequency
* Remove OTUs with relatively high paralogy frequency
* Mask monophylies by keepipng the longest sequence or the sequence with the shortest pairwise distance
* Exclude individual OTUs entirely
* Root trees using outgroup or midpoint rooting
* Get rid of OTUs with sequences that display relatively high pairwise distance
* Measure the impact of individual OTUs using taxon jackknifing

## Installation

This software runs under both Python 3 and 2.7. There are no external
dependencies, but [the plotting library Matplotlib](https://matplotlib.org/)
can be installed for generating paralog frequency plots.

You can install PhyloPyPruner using pip.

```bash
pip install --user phylopypruner
```

## Usage

To get a list of options, run the software without any arguments or use the
`-h` or `--help` flag. PhyloPyPruner requires either a corresponding multiple
sequence alignment (MSA) in FASTA format and a Newick tree or, the path to a
directory containing multiple trees and alignments.

**Example 1.** Providing a single corresponding tree and alignment. In this
case monophyletic masking will be performed by choosing the sequence with the
shorter pairwise distance to its sister group and paralogy pruning will be done
using the largest subtree (LS) algorithm.

```bash
python -m phylopypruner --msa <filename>.fas --tree <filename>.tre
```

**Example 2.** Run PhyloPyPruner for every MSA and tree pair within the
directory in `<path>`. Don't include orthologs with fewer than 10 OTUs, remove
sequence shorter than 100 positions, collapse nodes with a support value lower
than 80% into polytomies, remove branches that are 5 times longer than the
standard deviation of all branch lengths and remove OTUs with a paralogy
frequency that is larger than 5 times the standard deviation of the paralogy
frequency for all OTUs.

```bash
python -m phylopypruner --dir <path> --min-taxa 10 --min-len 100 --min-support
80 --trim-lb 5 --trim-freq-paralogs 5
```

**Example 3.** Run PhyloPyPruner for every MSA and tree pair within the
directory in `<path>`. Mask monophylies by choosing the longest sequence, prune
paralogs using the maximum inclusion (MI) algorithm, remove OTUs with sequences
with an average pairwise distance that is 10 times larger than the standard
deviation of the average pairwise distance of the sequences for all OTUs and
generate statistics for the removal of OTUs using taxon jackknifing.

```bash
python -m phylopypruner --dir <path> --mask longest --prune MI --trim-divergent
10 --jackknife
```

>>>
**Note:** Taxon jackknifing multiplies the execution time by the amount of OTUs
available within each input alignment.
>>>

FASTA descriptions and Newick names must match and has to be in one of the
following formats: `OTU|ID` or `OTU@ID`, where `OTU` is the operational
taxonomical unit (usually the species) and `ID` is a unique annotation or
sequence identifier. For example: `>Meiomenia_swedmarki|Contig00001_Hsp90`.
Sequence descriptions and tree names are not allowed to deviate from each
other. Sequence data needs to be [valid IUPAC nucleotide or amino acid
sequences](https://www.bioinformatics.org/sms/iupac.html).

## Output files

The following files are generated after running this program.

```
<output directory>/
├── <timestamp>_ppp_summary.csv
├── <timestamp>_ppp_ortho_stats.csv
├── <timestamp>_ppp_run.log
├── <timestamp>_ppp_paralog_freq.csv
├── <timestamp>_ppp_paralog_freq.png*
└── <timestamp>_orthologs/
│   ├── 1_pruned.fas
│   ├── 2_pruned.fas
│   ├── 3_pruned.fas
│   └── 4_pruned.fas
...
```

If `<output directory>` has not been specified by the `--output` flag, then
output files will be stored within the same directory as the input alignment
file(s). See the [Output files
section](https://gitlab.com/fethalen/phylopypruner/wikis/Output-Files) within
[the Wiki](https://gitlab.com/fethalen/phylopypruner/wikis/home) for a more
detailed
[explanation](https://gitlab.com/fethalen/phylopypruner/wikis/Output-Files#explanation)
of each individual output file.

\* – only produced if [Matplotlib](https://matplotlib.org/) is installed

© [Kocot Lab](https://www.kocotlab.com/) 2018
