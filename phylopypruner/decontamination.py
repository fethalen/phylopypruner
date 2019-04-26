"Module for dealing with contamination-like issues."

from __future__ import print_function
from __future__ import absolute_import
import sys
import copy
import datetime
from collections import defaultdict
from itertools import combinations
from functools import partial
from multiprocessing import Pool
from multiprocessing import cpu_count
from phylopypruner import filtering
from phylopypruner.summary import Summary
from phylopypruner.prune_paralogs import prune_paralogs

TIMESTAMP = datetime.datetime.now().strftime("%Y-%m-%d")
SUBCLADE_STATS_FILE = "/subclade_stats.csv"

def _exclude_and_rerun(taxon, summary, pruning_method, min_taxa, outgroup, dir_out):
    summary_copy = copy.deepcopy(summary)

    resample_summary = Summary()
    for log in summary_copy.logs:
        log_resampled = _resample(log, taxon, pruning_method, min_taxa,
                                    outgroup, dir_out)
        if log_resampled:
            resample_summary.logs.append(log_resampled)

    return resample_summary

def _resample(log, excluded, pruning_method, min_taxa, outgroup, dir_out):
    resample_log = copy.deepcopy(log)
    tree_excluded = filtering.exclude(resample_log.masked_tree, excluded)
    resample_log.msas_out = []
    if not tree_excluded:
        return None
    resample_log.settings.exclude = excluded
    resample_log.orthologs = prune_paralogs(pruning_method,
                                            tree_excluded,
                                            min_taxa,
                                            outgroup)
    resample_log.get_msas_out(dir_out)
    return resample_log

def jackknife(summary, dir_out, threads):
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
    taxa = summary.otus()
    total = len(taxa)
    resamples = set()
    # reuse the settings from the first log in the summary
    log = summary.logs[0]
    pruning_method = log.settings.prune
    min_taxa = log.settings.min_taxa
    outgroup = log.settings.outgroup
    pool = Pool(threads)

    part_jackknife = partial(
        _exclude_and_rerun, summary=copy.deepcopy(summary),
        pruning_method=pruning_method, min_taxa=min_taxa, outgroup=outgroup,
        dir_out=dir_out)

    for index, resample_summary in enumerate(
            pool.imap_unordered(part_jackknife, taxa), 1):
        print("{}==>{} jackknife resampling ({}/{} subsamples)".format(
            "\033[34m", "\033[0m", index, total), end="\r")
        sys.stdout.flush()
        resamples.add(resample_summary)
    pool.terminate()
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

def _rerun_wo_otu(log, otus, dir_out):
    log_copy = copy.deepcopy(log)
    pruning_method = log_copy.settings.prune
    min_taxa = log_copy.settings.min_taxa
    outgroup = log_copy.settings.outgroup
    tree = log_copy.masked_tree
    tree_excluded = filtering.exclude(tree, list(otus))
    log_copy.msas_out = []
    log_copy.settings.exclude = list(otus)
    if not tree_excluded:
        return None
    log_copy.orthologs = prune_paralogs(pruning_method,
                                        tree_excluded,
                                        min_taxa,
                                        outgroup)
    log_copy.get_msas_out(dir_out)
    return log_copy

def exclude_otus(summary, otus):
    """Exclude the provided OTUs from the provided Summary object.

    Parameters
    ----------
    summary : Summary object
        Prune MSAs from this summary.
    otus : list
        Remove the OTUs within this list from the Summary object.

    Returns
    -------
    summary : Summary object
        Input summary with the provided OTUs excluded.
    """
    for log in summary.logs:
        for msa in log.msas_out:
            for sequence in msa.sequences:
                if sequence.otu in otus:
                    msa.sequences.remove(sequence)

    return summary

def exclude_genes(summary, msas):
    """Exclude the multiple sequence alignments (MSAs) within the provided list
    from the provided Summary object.

    Parameters
    ----------
    summary : Summary object
        Prune MSAs from this summary.
    msas : list
        Remove the MultipleSequenceAlignment objects within this list from the
        Summary object.

    Returns
    -------
    summary : Summary object
        Input summary with the provided MSAs excluded.
    """
    for log in summary.logs:
        for msa in log.msas_out:
            if msa in msas:
                log.msas_out.remove(msa)

    return summary

def prune_by_exclusion(summary, otus, dir_out, threads, homolog_stats):
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
    excluded_str = "output_{}".format("+".join(otu for otu in otus))
    excluded_str += "_excluded"

    part_rerun = partial(_rerun_wo_otu, otus=otus, dir_out=dir_out)
    pool = Pool(threads)

    for index, log_copy in enumerate(
            pool.imap_unordered(part_rerun, summary_copy.logs), 1):
    # for index, log in enumerate(summary_copy.logs, 1):
        print("{}==>{} paralogy pruning with OTUs removed ({}/{} trees)".format(
            "\033[34m", "\033[0m", index, alignments_count), end="\r")
        sys.stdout.flush()

        if log_copy:
            summary_out.logs.append(log_copy)

    pool.terminate()
    print("")
    report = summary_out.report(excluded_str, dir_out, homolog_stats)

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
    threshold = _std(list(paralog_freq.values())) * factor
    otus_above_threshold = list()

    for otu in paralog_freq:
        if paralog_freq[otu] > threshold:
            otus_above_threshold.append(otu)

    if not otus_above_threshold:
        print("OTUs with high paralogy frequency: none")
        return otus_above_threshold

    print("OTUs with high paralogy frequency: " +
          ", ".join(otu for otu in otus_above_threshold))

    return otus_above_threshold

def trim_divergent(node, divergence_threshold=0.25, include=[]):
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

    if include:
        for otu in include:
            if otu in nodes_to_remove:
                nodes_to_remove.remove(otu)

    node.remove_nodes(nodes_to_remove)

    return otus_above_threshold

def score_monophyly(summary, taxonomic_groups, dir_out):
    """Takes a Summary object, a list of TaxonomicGroup objects, and the path
    to the output directory as an input. For each ortholog (output alignment)
    within the Summary object, analyse the monophyly of each group defined
    within the list of TaxonomicGroup objects by counting how many times each
    group forms a monophyletic group, which OTUs are present and which OTUs are
    "invading" a group. Output a CSV file where one axis is the OTUs and
    the other is the various groups.

    Parameters
    ----------
    summary : Summary object
        The orthologs (output alignments) from this Summary object are used in
        this analysis.
    taxonomic_groups : list
        This objects contains the different groups: their name and the OTUs
        within the group. Each item in the list should be a TaxonomicGroup
        object.
    dir_out : str
        The directory of where the output statistics file is written to.
    """
    otus = set()
    otu_scores_ordered = list()
    group_scores = defaultdict(int)
    row = defaultdict(int)
    group_names = set([group.name for group in taxonomic_groups])
    trees = set()

    print("{}==>{} calculating monophyly scores".format("\033[34m", "\033[0m"))

    # compile sets of trees and OTUs
    for log in summary.logs:
        for ortholog in log.orthologs:
            trees.add(ortholog)
            for otu in ortholog.iter_otus():
                otus.add(otu)

    otus = list(otus)

    for group in taxonomic_groups:
        group_scores[group.name] = [0 for otu in otus]

    otu_scores_ordered = [0 for otu in otus]

    for tree in trees:
        for group in taxonomic_groups:
            ingroups = tree.outgroups_present(group.otus)
            most_inclusive_branch = None
            most_ingroups = 0
            no_of_outgroups = 0

            # only consider cases where 2, or more, members are present
            if len(ingroups) <= 1:
                continue

            for branch in tree.iter_branches():
                otus_in_branch = set(branch.iter_otus())
                ingroups_in_branch = otus_in_branch.intersection(ingroups)
                outgroups_in_branch = otus_in_branch.difference(ingroups_in_branch)

                if len(ingroups_in_branch) < 2 or len(outgroups_in_branch) > 2:
                    continue

                # Find the branch that maximizes the amount of ingroup OTUs and
                # minimizes the amount of outgroup OTUs.
                if len(ingroups_in_branch) > most_ingroups or \
                    (len(ingroups_in_branch) == most_ingroups and \
                    len(outgroups_in_branch) < no_of_outgroups):
                    most_inclusive_branch = branch
                    most_ingroups = len(ingroups_in_branch)
                    no_of_outgroups = len(outgroups_in_branch)

            if most_inclusive_branch:
                for otu in most_inclusive_branch.iter_otus():
                    # Add +1 for each time OTU forms a monophyletic group with
                    # group X.
                    index = otus.index(otu)
                    group_scores[group.name][index] += 1

                    if otu in ingroups:
                        # OTU present and forms a monophyletic group
                        otu_scores_ordered[index] += 2
                        ingroups.remove(otu)
                        # group_scores[otus.index(otu)] += 1
                    else:
                        # OTU is present but invades another group
                        otu_scores_ordered[index] -= 2
                # Score the OTUs which were present but did not form a
                # monophyletic group.
                for otu in ingroups:
                    otu_scores_ordered[index] -= 2

    # Order the score so that there is one OTU for each row, instead of for
    # each column.
    scores_otu_row = []
    for index, score in enumerate(otu_scores_ordered):
        scores_otu_row.append([otus[index], score])

    for group in group_scores:
        for index, score in enumerate(group_scores[group]):
            scores_otu_row[index].append(group_scores[group][index])

    # write the results to a CSV file.
    with open(dir_out + SUBCLADE_STATS_FILE, "w") as subclades_file:
        subclades_file.write("otu;monophylyScore;" + ";".join([group for group in group_names]) +
                "\n")
        for row in scores_otu_row:
            subclades_file.write(";".join([str(item) for item in row]) +
                    "\n")
