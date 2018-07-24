PhyloPyPruner
-------------

PhyloPyPruner is a phylogenetic tree-based orthology inference program, that
implements a variety of filters and algorithms, to select putative orthologs
from a set trees and sequences.

PhyloPyPruner succeeds
[PhyloTreePruner](https://sourceforge.net/projects/phylotreepruner/), in that
it implements a superset of the functionality found in mentioned software.

### Input data

The `--msa` flag, short for multiple sequence alignment (MSA), points to one or
more alignments, or a directory of alignments. The `-tree` flag behaves
similarly, in that it can take a set of trees, or a directory of trees.

```
./phylopypruner --msa test.fa --tree test.nw
```

### Algorithms and parameters

Orthology inference can be divided into several steps. Here, each step is
implemented as a stand-alone Python module, that are all accessible from a
single-user interface, for maximum flexibility.

#### Homology inference

`infer_homologs.py` module.

#### Filtering-steps

| **Command**                   | **Meaning**                                                                           |
|---------------------------|-----------------------------------------------------------------------------------|
| `--minimum-positions cutoff` | Trim sequences shorter than `cutoff` positions                                    |
| `--minimum-taxa cutoff`  | Remove multiple sequence alignments where less than `cutoff` taxa are represented |

#### Monophyly masking/isoform pruning

#### Alignment

Align using [MAFFT](https://mafft.cbrc.jp/alignment/software/) and trim
alignment using
[Gblocks](http://molevol.cmima.csic.es/castresana/Gblocks.html).

#### Tree Inference

#### Rooting

`rooting.py` module. No rooting is performed, if no rooting method is
specified.

The root represents the last common ancestor (LCA) from which all other
operational taxonomic units (OTUs) in the tree are descended.

**Outgroup**

One or more OTUs – that is, outgroup(s) – are selected based on prior
phylogenetic knowledge about the OTUs under study. The branch that harbors all
other OTUs – that is, ingroups – is assumed to be the root. A wrong assumption
about the relationship between the outgroups and ingroups will lead to a wrong
tree topology.

**Midpoint**

The distance between all OTU pairs is calculated by adding the length of the
intervening branches together and the root is placed in the middle of the
longest path.

**Relaxed molecular-clock**

**Minimal ancestor deviation (MAD)**

#### Paralog pruning

`prune_paralogs.py` module. Recursive maximum inclusion (RMI) is used, if no
algorithm is specified.

| **Command**                   | **Meaning**                                                                           |
|---------------------------|-----------------------------------------------------------------------------------|
| `--collapse-nodes support` | Collapse nodes with a support value below
`support`                                    |

**Maximum inclusion (MI)**

**Iterative maximum inclusion (IMI)**

**Monophyletic outgroups (MO)**

**Rooted tree (RT)**

**One-to-one orthologs (1to1)**

#### Contamination screening

© [Kocot Lab](https://www.kocotlab.com/) 2018
