"Needed for multiprocessing (see: https://bugs.python.org/issue25053)"

import copy
import sys
from . import decontamination
from . import prune_paralogs
from . import mask_monophylies
from . import filtering
from . import fasta
from . import log
from . import newick
from . import report
from . import root

version = "1.1.2"

def validate_input(msa, tree, tree_path):
    "Test to see if MSA and tree entries matches."
    descriptions = list(msa.iter_descriptions())
    names = list(tree.iter_names())

    if set(descriptions).intersection(names) < set(descriptions):
        print("example tree names:", names[:2], file=sys.stderr)
        print("example sequences:", descriptions[:2], file=sys.stderr)
        report.error("MSA names don't match tree \n   {}\n   {}".format(
            msa.filename, tree_path))


def run(settings, msa, tree):
    """
    Takes a dictionary that specifies a set of settings as an input. The
    dictionary should contain the following keys: msa, tree, min_taxa, min_len,
    min_support, trim_lb, outgroup, root, mask, and prune. Runs PhyloPyPruner
    once using these settings.
    """
    # test to see if the entries in the MSA and tree matches
    validate_input(msa, tree, settings.nw_file)

    run_log = log.Log(version, msa, tree, settings)

    # exclude taxa within the list settings.exclude
    if settings.exclude:
        tree = filtering.exclude(tree, settings.exclude)

    if filtering.too_few_otus(tree, settings.min_taxa):
        return run_log

    # remove short sequences
    if settings.min_len:
        run_log.trimmed_seqs = filtering.trim_short_seqs(
            msa, tree, settings.min_len)

    if filtering.too_few_otus(tree, settings.min_taxa):
        return run_log

    # trim long branches
    if settings.trim_lb:
        run_log.lbs_removed = list(filtering.prune_long_branches(
            tree, settings.trim_lb))

    if filtering.too_few_otus(tree, settings.min_taxa):
        return run_log

    # trim zero length branches
    if settings.trim_zero_len:
        run_log.ultra_short_branches = filtering.trim_zero_len_branches(
            tree, settings.trim_zero_len)

    if filtering.too_few_otus(tree, settings.min_taxa):
        return run_log

    # collapse weakly supported nodes into polytomies
    if settings.min_support:
        run_log.collapsed_nodes = filtering.collapse_nodes(
            tree, settings.min_support)

    # mask monophyletic groups
    if settings.mask:
        if settings.mask == "pdist":
            tree, masked_seqs = mask_monophylies.pairwise_distance(tree)
        elif settings.mask == "longest":
            tree, masked_seqs = mask_monophylies.longest_isoform(msa, tree)
        run_log.monophylies_masked = masked_seqs

    tree, masked_seqs = mask_monophylies.pairwise_distance(tree)

    # trim divergent sequences
    if settings.trim_divergent:
        run_log.divergent, run_log.divergent_removed = decontamination.trim_divergent(
            tree, settings.trim_divergent, settings.include)

    if filtering.too_few_otus(tree, settings.min_taxa):
        return run_log

    # root by outgroup
    rooted = False   # True if outgroup rooting is successful
    if settings.outgroup:
        if not settings.prune == "MO":
            tree, rooted = root.outgroup(tree, settings.outgroup)

    # root the tree by midpoint rooting
    if not rooted and settings.root:
        if settings.root == "midpoint":
            tree = root.midpoint(tree)

    # mask monophyletic groups
    if settings.mask:
        if settings.mask == "pdist":
            tree, masked_seqs = mask_monophylies.pairwise_distance(tree)
        elif settings.mask == "longest":
            tree, masked_seqs = mask_monophylies.longest_isoform(msa, tree)
        run_log.monophylies_masked.update(masked_seqs)

    run_log.masked_tree = copy.deepcopy(tree)

    # get a list of paralogs
    run_log.paralogs = tree.paralogs()

    # prune paralogs
    run_log.orthologs = prune_paralogs.prune_paralogs(settings.prune, tree,
                                   settings.min_taxa, settings.outgroup)

    if settings.force_inclusion:
        run_log.orthologs = filtering.force_inclusion(
            run_log.orthologs, settings.force_inclusion)

    return run_log


def run_for_file_pairs(corr_files, settings, dir_in, dir_out):
    settings.fasta_file, settings.nw_file = corr_files
    return get_orthologs(settings, dir_in, dir_out)


def get_orthologs(settings, directory="", dir_out=None):
    fasta_path = "{}{}".format(directory, settings.fasta_file)
    nw_path = "{}{}".format(directory, settings.nw_file)
    msa = fasta.read(fasta_path)
    nw_file = newick.read(nw_path)
    log = run(settings, msa, nw_file)
    log.get_msas_out(dir_out)
    log.report(dir_out)
    return log
