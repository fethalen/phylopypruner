"Module for running various decontamination tools."

from __future__ import print_function
import sys
import copy
import datetime
from collections import defaultdict
from itertools import combinations
import filtering
from summary import Summary
from prune_paralogs import prune_paralogs

TIMESTAMP = datetime.datetime.now().strftime("%Y-%m-%d")

def jackknife(summary, dir_out):
    resampling = set(combinations(summary.otus(), len(summary.otus()) - 1))
    total = len(resampling)
    resamples = set()

    for index, combination in enumerate(resampling, 1):
        summary_copy = copy.deepcopy(summary)
        resample_summary = Summary()
        excluded = list(summary.otus().difference(combination))

        sys.stdout.flush()
        print("jackknife resampling ({}/{} subsamples)".format(
            index, total), end="\r")

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
    """ Return the sample arithmetic mean of data, a sequence of real-valued
    numbers. The average of the empty list, '[]', is 0.
    """
    return float(sum(data)) / max(len(data), 1)

def _sdm(data):
    """ Return the sum of square deviations of data.
    """
    return sum((x - _mean(data))**2 for x in data)

def _std(data):
    "Return the population standard deviation of data."
    if len(data) < 2:
        raise ValueError('variance requires at least two data points')
    return (_sdm(data) / len(data)) ** 0.5

def prune_by_exclusion(summary, otus, dir_out):
    """
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

def trim_freq_paralogs(summary, factor, paralog_freq):
    treshold = _std(list(paralog_freq.values())) * factor
    rogue_taxa = set()

    for otu in paralog_freq:
        if paralog_freq[otu] > treshold:
            rogue_taxa.add(otu)

    rogue_taxa_str = "+".join(otu for otu in rogue_taxa)
    rogue_taxa_str += "_excluded"

    if not rogue_taxa:
        print("OTUs with frequent paralogs: none")
        return

    print("OTUs with frequent paralogs: " +
          ", ".join(otu for otu in rogue_taxa))

    return rogue_taxa

def trim_divergent_otus(summary, factor):
    """
    Takes a Summary object, an integer and the path to an output directory as
    an input. Calculates the average pairwise distance for each OTU within
    <summary>. OTUs with an average pairwise distance that is <factor> times
    the standard deviation of all average pairwise distances are removed from
    the trees within the Summary object's logs. All trees are then pruned using
    the same method that was provided in the summary.
    """
    otu_dists = defaultdict(float) # OTU is key, distance is value
    dists_count = defaultdict(int) # OTU is key, frequency is value
    otus_above_threshold = list()
    total = len(summary)

    for index, log in enumerate(summary.logs, 1):

        sys.stdout.flush()
        print("calculating pairwise distances ({}/{} trees)".format(
            index, total), end="\r")

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
