What is PhyloPyPruner?
----------------------

PhyloPyPruner is a Python package for refining the output of a graph-based
orthology inference approach such as
[OrthoMCL](https://www.ncbi.nlm.nih.gov/pubmed/12952885),
[OrthoFinder](https://www.ncbi.nlm.nih.gov/pubmed/26243257) or
[HaMStR](https://www.ncbi.nlm.nih.gov/pubmed/19586527). Similar to other
tree-based orthology inference methods (e.g.,
[PhyloTreePruner](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3825643/),
[UPhO](https://academic.oup.com/mbe/article/33/8/2117/2578877),
[Agalma](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3840672/) and
[Phylogenomic Dataset
Reconstruction](https://www.ncbi.nlm.nih.gov/pubmed/25158799)), it
uses phylogenetic trees in order to obtain genes that are 1:1 orthologous. In
addition to filters and algorithms seen in pre-existing tools, this package
provides new methods for differentiating [contamination-like
sequences](https://gitlab.com/fethalen/phylopypruner/wikis/About-PhyloPyPruner#contamination-like-issues-)
from paralogs.

PhyloPyPruner is under active development and I would appreciate it if you try
this software on your own data and [leave
feedback](mailto:felix.thalen@uni-goettingen.de).

![proteomes2orthologs](https://gitlab.com/fethalen/phylopypruner/raw/master/doc/images/proteomes2orthologs.png)

**Figure 1.** An overview of a tree-based orthology inference approach. First,
sequences from different species are clustered together based on an
all-versus-all BLAST, followed by Markov clustering. Each node in the cluster
corresponds to a sequence and each edge corresponds to a similarity score.
Homologous groups are then aligned and a phylogenetic tree is inferred from the
alignment. From this tree, orthologous groups can be identified and paralogs
are pruned away.

## Quick installation

The easiest way to install PhyloPyPruner is by using the package manager for
Python, [pip](https://pypi.org/project/pip/):

```bash
pip install phylopypruner # install for all users
pip install --user phylopypruner # install for the current user only
```

Once installed, the program is located within `$HOME/.local/bin`. Depending on
your OS, you might have to add the directory to your `$PATH` to avoid typing
the entire path. Once in your path, you run the program like this:

```bash
phylopypruner
```

## [Documentation](https://gitlab.com/fethalen/phylopypruner/wikis)

* [About PhyloPyPruner](https://gitlab.com/fethalen/phylopypruner/wikis/about-phylopypruner)
* Tutorial
* [Installation](https://gitlab.com/fethalen/phylopypruner/wikis/installation)
* [Input data](https://gitlab.com/fethalen/phylopypruner/wikis/input-data)
* [Output files](https://gitlab.com/fethalen/phylopypruner/wikis/output-files)
* [Methods](https://gitlab.com/fethalen/phylopypruner/wikis/methods)
* [Options](https://gitlab.com/fethalen/phylopypruner/wikis/options)

## Cite

Our manuscript is still in preparation, it will be posted here once a preprint
of the article is available.

© [Kocot Lab](https://www.kocotlab.com/) 2018
