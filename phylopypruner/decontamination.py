"Module for running various decontamination tools."

from __future__ import print_function
import sys
import copy
import datetime
from itertools import combinations
import filtering
from summary import Summary
from prune_paralogs import prune_paralogs

TIMESTAMP = datetime.datetime.now().strftime("%Y-%m-%d")

def otus(summary):
    "Returns a set of all OTUs within all of the trees in the Summary object."
    otus_in_summary = set()
    for log in summary.logs:
        otus = set(list(log.masked_tree.iter_otus()))
        if otus:
            otus_in_summary.update(otus)

    return otus_in_summary

def jackknife(summary, dir_out):
    otus_in_summary = otus(summary)
    otu_count = len(otus_in_summary)

    resampling = set(combinations(otus_in_summary, otu_count - 1))
    total = len(resampling)
    resamples = set()

    for index, combination in enumerate(resampling, 1):
        summary_copy = copy.deepcopy(summary)
        resample_summary = Summary()
        excluded = list(otus_in_summary.difference(combination))

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

def trim_freq_paralogs(summary, factor, paralog_freq, dir_out):
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


    summary_copy = copy.deepcopy(summary)
    rtr_summary = Summary()
    total = len(summary_copy.logs)

    for index, log in enumerate(summary_copy.logs, 1):
        sys.stdout.flush()
        print("paralogy pruning with frequent paralogs removed ({}/{} file \
pairs)".format(index, total), end="\r")

        rtr_log = copy.deepcopy(log)
        pruning_method = rtr_log.settings.prune
        min_taxa = rtr_log.settings.min_taxa
        outgroup = rtr_log.settings.outgroup
        rtr_tree = rtr_log.masked_tree
        tree_excluded = filtering.exclude(rtr_tree, list(rogue_taxa))
        rtr_log.msas_out = []
        rtr_log.settings.exclude = list(rogue_taxa)
        if not tree_excluded:
            continue
        rtr_log.orthologs = prune_paralogs(pruning_method,
                                           tree_excluded,
                                           min_taxa,
                                           outgroup)
        rtr_log.get_msas_out(dir_out)
        rtr_summary.logs.append(rtr_log)

    print("")
    rtr_summary_report = rtr_summary.report(rogue_taxa_str, dir_out)

    return rtr_summary_report, rtr_summary
