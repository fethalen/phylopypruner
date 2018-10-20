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

def trim_divergent(node, divergence_threshold=0.25):
    """For each OTU with more than one sequence present in the provided node:
    calculate the ratio of the maximum pairwise distance of the sequences
    within the OTU compared to the average pairwise distance for that OTU
    compared to every other sequence. Delete every sequence from that OTU
    entirely if the ratio exceeds the divergence threshold.

    Parameters
    ----------
    node : TreeNode object
        The node that you want wish to delete divergent sequences from.
    divergence_threshold : float
        Divergence threshold in percent.

    Returns
    ------
    seqs_above_threshold : list
        List of OTUs above the established threshold.
    """
    # maximum pairwise distance within OTUs ({OTU: max_pdist})
    in_otus_max_dist = defaultdict(float)
    # pairwises pairwise distance between one OTU's sequences and other OTU's
    # sequences ({OTU: [dist_1, dist_2, ...]})
    out_otus_dists = defaultdict(list)
    # average pairwise distance between one OTU's sequences and another OTU's
    # sequences ({OTU: avg_pdist})
    out_otus_avg_dist = dict()
    # ratio between the maximum pairwise distance of the in-OTUs and the
    # average pairwise distance of the out-OTUs
    in_out_ratio = dict()
    otus_above_threshold = list()
    nodes_to_remove = set()
    otus_removed = 0

    for paralog in node.paralogs():
        for leaf in node.iter_leaves():
            if paralog is leaf:
                continue

            if paralog.otu() == leaf.otu():
                if (not in_otus_max_dist[paralog.otu()] or
                        paralog.distance_to(leaf) >
                        in_otus_max_dist[paralog.otu()]):
                    in_otus_max_dist[paralog.otu()] = paralog.distance_to(leaf)
            else:
                out_otus_dists[paralog.otu()].append(paralog.distance_to(leaf))

    for otu in out_otus_dists:
        out_otus_avg_dist[otu] = sum(out_otus_dists[otu]) / float(len(out_otus_dists[otu]))

    for otu in out_otus_avg_dist:
        in_out_ratio = in_otus_max_dist[otu] / out_otus_avg_dist[otu]
        if in_out_ratio > divergence_threshold:
            otus_above_threshold.append(otu)

    for leaf in node.iter_leaves():
        if leaf.otu() in otus_above_threshold:
            nodes_to_remove.add(leaf)
            otus_removed += 1

    node.remove_nodes(nodes_to_remove)

    return otus_removed
