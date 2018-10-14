"Module for dealing with contamination-like issues."

from __future__ import print_function
from __future__ import absolute_import
import sys
import copy
import datetime
from collections import defaultdict
from itertools import combinations
from phylopypruner import filtering
from phylopypruner.summary import Summary
from phylopypruner.prune_paralogs import prune_paralogs

TIMESTAMP = datetime.datetime.now().strftime("%Y-%m-%d")

def jackknife(summary, dir_out):
    """Exclude each OTUs within the summary, one by one, perform paralogy
    pruning and output summary statistics of the output alignments for each
    subsample.

    Parameters
    ----------
    summary : Summary object
        Perform paralogy pruning on the Log object's masked tree attribute, for
        each Log object within this Summary's logs attribute.
    dir_out : str
        Write the statistics for each case to the summary file within this
        directory.

    Returns
    -------
    None
    """
    resampling = set(combinations(summary.otus(), len(summary.otus()) - 1))
    total = len(resampling)
    resamples = set()

    for index, combination in enumerate(resampling, 1):
        summary_copy = copy.deepcopy(summary)
        resample_summary = Summary()
        excluded = list(summary.otus().difference(combination))

        print("jackknife resampling ({}/{} subsamples)".format(
            index, total), end="\r")
        sys.stdout.flush()

        for log in summary_copy.logs:
            resample_log = copy.deepcopy(log)
            pruning_method = log.settings.prune
            min_taxa = log.settings.min_taxa
            outgroup = log.settings.outgroup
            resample_tree = copy.copy(log.masked_tree)
            tree_excluded = filtering.exclude(resample_tree, excluded)
            resample_log.msas_out = []
            if not tree_excluded:
                continue
            resample_log.settings.exclude = excluded
            resample_log.orthologs = prune_paralogs(pruning_method,
                                                    tree_excluded,
                                                    min_taxa,
                                                    outgroup)
            resample_log.get_msas_out(dir_out)
            resample_summary.logs.append(resample_log)

        resamples.add(resample_summary)
    print("")

    for summary in resamples:
        excluded = summary.logs[0].settings.exclude[0]
        summary.report("{}_excluded".format(excluded), dir_out)

def _mean(data):
    """Returns the sample arithmetic mean of data. 0 is returned if an empty
    list was provided.

    Parameters
    ----------
    data : list of floats

    Returns
    _______
    out: float
        The sample arithmetic mean of data.
    """
    return float(sum(data)) / max(len(data), 1)

def _sdm(data):
    """Returns the squared deviations from the mean (SDM) of data.

    Parameters
    ----------
    data : list of floats

    Returns
    -------
    out : float
        The sum of square deviations of data.
    """
    return sum((x - _mean(data))**2 for x in data)

def _std(data):
    """Return the population standard deviation of data.

    Parameters
    ----------
    data : list of floats

    Returns
    -------
    out : float
        The population standard deviation of data.
    """
    if len(data) < 2:
        raise ValueError('variance requires at least two data points')
    return (_sdm(data) / len(data)) ** 0.5

def prune_by_exclusion(summary, otus, dir_out):
    """Exclude the OTUs within the provided list of OTUs from the masked trees
    within summary, perform paralogy pruning and output statistics and
    alignments for each ortholog recovered.

    Parameters
    ----------
    summary : Summary object
        Prune masked trees within the logs of this summary.
    otus : list of strings
        Exclude OTUs within this list from the trees in the summary.
    dir_out : str
        Output statistics to this directory.

    Returns
    -------
    summary_out : Summary object
        A new summary that was generated after performing paralogy pruning on
        the masked trees within the input summary with the OTUs within the
        provided list removed.
    report : str
        Printable statistics of the summary_out

    Takes a Summary object, a list of OTUs and the path to the output directory
    as an input. Returns a new Summary object that is the summary after
    paralogy pruning with the OTUs within the list excluded.
    """
    # creating a copy of the summary prevents making changes to the trees in
    # that summary
    summary_copy = copy.deepcopy(summary)
    summary_out = Summary()
    alignments_count = len(summary_copy.logs)
    excluded_str = "+".join(otu for otu in otus)
    excluded_str += "_excluded"

    for index, log in enumerate(summary_copy.logs, 1):
        sys.stdout.flush()
        print("paralogy pruning with OTUs removed ({}/{} trees)".format(
            index, alignments_count), end="\r")

        log_copy = copy.deepcopy(log)
        pruning_method = log_copy.settings.prune
        min_taxa = log_copy.settings.min_taxa
        outgroup = log_copy.settings.outgroup
        tree = log_copy.masked_tree
        tree_excluded = filtering.exclude(tree, list(otus))
        log_copy.msas_out = []
        log_copy.settings.exclude = list(otus)
        if not tree_excluded:
            continue
        log_copy.orthologs = prune_paralogs(pruning_method,
                                            tree_excluded,
                                            min_taxa,
                                            outgroup)
        log_copy.get_msas_out(dir_out)
        summary_out.logs.append(log_copy)
    print("")

    report = summary_out.report(excluded_str, dir_out)
    return summary_out, report

def trim_freq_paralogs(factor, paralog_freq):
    """Returns a set of OTUs with a paralogy frequency that is factor times
    larger than the standard deviation of the paralogy frequency of all OTUs.

    Parameters
    ----------
    factor : float
        Set the threshold to be this float multiplied by the standard deviation
        of the paralogy frequency of all OTUs.
    paralog_freq : dictionary
        Paralogy frequency for each OTU, where key is OTU and paralogy
        frequency is value.

    Returns
    _______
    otus_above_threshold : list
        A set of OTUs with a paralogy frequency above the threshold.
    """
    treshold = _std(list(paralog_freq.values())) * factor
    otus_above_threshold = list()

    for otu in paralog_freq:
        if paralog_freq[otu] > treshold:
            otus_above_threshold.append(otu)

    if not otus_above_threshold:
        print("OTUs with frequent paralogs: none")
        return

    print("OTUs with frequent paralogs: " +
          ", ".join(otu for otu in otus_above_threshold))

    return otus_above_threshold

def trim_divergent_otus(summary, factor):
    """Return a list of OTUs in the summary with an average pairwise distance
    that is larger than factor times the standard deviation of the pairwise
    distance of all OTUs.

    Parameters
    ----------
    summary : Summary object
        The pairwise distance is derived from the trees within this Summary
        object.
    factor : float
        Set the threshold to be this float times the standard deviation of the
        pairwise distance for all OTUs.

    Returns
    ------
    otus_above_threshold : list
        List of OTUs above the established threshold.
    """
    otu_dists = defaultdict(float) # OTU is key, distance is value
    dists_count = defaultdict(int) # OTU is key, frequency is value
    otus_above_threshold = list()
    total = len(summary)

    for index, log in enumerate(summary.logs, 1):

        print("comparing pairwise distances ({}/{} trees)".format(
            index, total), end="\r")
        sys.stdout.flush()

        if not log.masked_tree:
            # case were a 'masked tree' was not created due to too few OTUs
            continue

        dists_in_tree = log.masked_tree.distances()
        for pair in dists_in_tree:
            leaf_a, leaf_b = pair
            otu_a = leaf_a.otu()
            otu_b = leaf_b.otu()
            otu_dists[otu_a] += dists_in_tree[pair]
            otu_dists[otu_b] += dists_in_tree[pair]
            dists_count[otu_a] += 1
            dists_count[otu_b] += 1
    print("")

    for otu in dists_count:
        otu_dists[otu] = round(otu_dists[otu] / dists_count[otu], 2)

    threshold = _std(list(otu_dists.values())) * factor

    for otu in otu_dists:
        if otu_dists[otu] > threshold:
            otus_above_threshold.append(otu)

    if otus_above_threshold:
        print("OTUs with an average pairwise distance above threshold: \
{}".format(", ".join(otus_above_threshold)))
    else:
        print("no OTUs with an average pairwise distance above threshold")

    return otus_above_threshold

def trim_divergent_seqs(node, factor):
    """Iterate over all multiple sequence alignments within the provided
    summary and remove sequences with a pairwise distance that is larger than
    factor times the standard deviation of the pairwise distance of all
    sequences within that single gene.

    Parameters
    ----------
    summary : Summary object
        The pairwise distance is derived from the trees within this Summary
        object.
    factor : float
        Set the threshold to be this float times the standard deviation of the
        pairwise distance for all OTUs.

    Returns
    ------
    seqs_above_threshold : list
        List of OTUs above the established threshold.
    """
    seqs_above_threshold = 0

    seq_dists = defaultdict(float) # leaf is key, distance is value
    seq_count = defaultdict(int) # leaf is key, presence is value
    avg_dists = defaultdict(float) # leaf is key, average distance is value
    nodes_to_remove = set()

    dists_in_tree = node.distances()
    for pair in dists_in_tree:
        leaf_a, leaf_b = pair
        seq_dists[leaf_a] += dists_in_tree[pair]
        seq_dists[leaf_b] += dists_in_tree[pair]
        seq_count[leaf_a] += 1
        seq_count[leaf_b] += 1

    for leaf in seq_dists:
        avg_dists[leaf] = round(seq_dists[leaf] / seq_count[leaf], 2)

    print(avg_dists.values())
    threshold = _std(list(avg_dists.values())) * factor

    for leaf in avg_dists:
        if avg_dists[leaf] > threshold:
            seqs_above_threshold += 1
            nodes_to_remove.add(leaf)

    node.remove_nodes(nodes_to_remove)

    return seqs_above_threshold
