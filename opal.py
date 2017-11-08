#!/usr/bin/env python

import os
import errno
import argparse
import os.path
from collections import defaultdict
import l1norm as l1
import binary_metrics as bm
import unifrac_distance as uf
import shannon as sh
import plots as pl
from utils import load_data
from utils import ProfilingTools as PF
from utils import constants as c
import pandas as pd
import numpy as np


def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def get_labels(labels, profiles_files):
    if labels:
        labels_list = [x.strip() for x in labels.split(',')]
        if len(labels_list) != len(profiles_files):
            exit('Number of labels does not match the number of files of profiles. Please check parameter -l, --labels.')
        return labels_list
    tool_id = []
    for profile_file in profiles_files:
        tool_id.append(profile_file.split('/')[-1])
    return tool_id


def compute_binary_metrics(query_profile, query_truth):
    all_metrics = bm.compute_tree_metrics(query_profile, query_truth)
    return all_metrics


def print_by_rank(output_dir, labels, pd_metrics):
    for rank in c.ALL_RANKS:
        # define ordering of rows, which is given my order of tool labels
        order_rows = labels
        # define ordering of columns, hard coded
        order_columns = [c.UNIFRAC, c.UNW_UNIFRAC, c.L1NORM, c.RECALL, c.PRECISION, c.TP, c.FP, c.FN, c.JACCARD, c.SHANNON_DIVERSITY, c.SHANNON_EQUIT]
        # subset to those information that either belong to the given rank or are rank independent, i.e. are unifrac values
        table = pd_metrics[(pd_metrics['rank'] == rank) | (pd_metrics['metric'].isin([c.UNIFRAC, c.UNW_UNIFRAC]))]
        # reformat the table with a pivot_table
        table = table.pivot_table(index='tool', columns='metric', values='value')
        # order table
        table = table.loc[order_rows, order_columns]
        # replace np.NaN with string "na" and write resulting table into a file
        table.fillna('na').to_csv(os.path.join(output_dir, rank + ".tsv"), sep='\t')


def print_by_tool(output_dir, pd_metrics):
    # define ordering of columns, hard coded
    order_columns = [c.UNIFRAC, c.UNW_UNIFRAC, c.L1NORM, c.RECALL, c.PRECISION, c.TP, c.FP, c.FN, c.JACCARD, c.SHANNON_DIVERSITY, c.SHANNON_EQUIT]
    for toolname, pd_metrics_tool in pd_metrics.groupby('tool'):
        table = pd_metrics_tool.pivot_table(index='rank', columns='metric', values='value')
        # little hack to carry unifrac over to every rank
        for unifrac_col in order_columns[:2]:
            table[unifrac_col] = pd_metrics_tool[pd_metrics_tool['metric'] == unifrac_col]['value'].values[0]
        # order table
        table = table.loc[c.ALL_RANKS, order_columns]
        # replace np.NaN with string "na" and write resulting table into a file
        table.fillna('na').to_csv(os.path.join(output_dir, toolname + ".tsv"), sep='\t')


def evaluate(gold_standard_file, profiles_files, labels, output_dir):
    shannon_list = []
    l1norm_list = []
    binary_metrics_list = []
    weighted_unifrac_list = []
    gs_rank_to_taxid_to_percentage = load_data.open_profile(gold_standard_file)
    gs_profile = PF.Profile(input_file_name=gold_standard_file)
    metrics = [c.UNIFRAC, c.L1NORM, c.RECALL, c.PRECISION, c.FP]

    for profile_file, label in zip(profiles_files, labels):
        rank_to_taxid_to_percentage = load_data.open_profile(profile_file)
        pf_profile = PF.Profile(input_file_name=profile_file)

        # Shannon
        shannon = sh.compute_shannon_index(rank_to_taxid_to_percentage)
        shannon_list.append(shannon)

        # L1 Norm
        l1norm = l1.compute_l1norm(gs_rank_to_taxid_to_percentage, rank_to_taxid_to_percentage)
        l1norm_list.append(l1norm)

        # Binary metrics
        binary_metrics = compute_binary_metrics(rank_to_taxid_to_percentage, gs_rank_to_taxid_to_percentage)
        binary_metrics_list.append(binary_metrics)

        # Unifrac
        unifrac = uf.compute_unifrac(gs_profile, pf_profile)
        weighted_unifrac_list.append(unifrac)

    rank_to_metric_to_toolvalues = defaultdict(dict)
    for rank in c.PHYLUM_SPECIES:
        for metric in metrics:
            rank_to_metric_to_toolvalues[rank][metric] = []
        max_fp = max([profile[rank].fp for profile in binary_metrics_list])
        for unifrac, l1norm, binary_metrics in zip(weighted_unifrac_list, l1norm_list, binary_metrics_list):
            rank_to_metric_to_toolvalues[rank][c.UNIFRAC].append(unifrac[0] / 16)
            rank_to_metric_to_toolvalues[rank][c.L1NORM].append(l1norm[rank] / 2.0 if rank in l1norm else None)
            rank_to_metric_to_toolvalues[rank][c.RECALL].append(binary_metrics[rank].recall if rank in binary_metrics else None)
            rank_to_metric_to_toolvalues[rank][c.PRECISION].append(binary_metrics[rank].precision if rank in binary_metrics else None)
            if max_fp > 0:
                rank_to_metric_to_toolvalues[rank][c.FP].append(binary_metrics[rank].fp / max_fp if rank in binary_metrics else None)
            else:
                rank_to_metric_to_toolvalues[rank][c.FP].append(binary_metrics[rank].fp if rank in binary_metrics else None)

    pl.plot(metrics,
         labels,
         rank_to_metric_to_toolvalues,
         output_dir,
         'spider_plot',
         ['b', 'g', 'r', 'k', 'm'])

    pl.plot([c.RECALL, c.PRECISION],
         labels,
         rank_to_metric_to_toolvalues,
         output_dir,
         'spider_plot_recall_precision',
         ['r', 'k'],
         grid_points=[0.2, 0.4, 0.6, 0.8, 1.0],
         fill=True)

    rank_to_shannon_gs = sh.compute_shannon_index(gs_rank_to_taxid_to_percentage)
    pl.plot_shannon(shannon_list, rank_to_shannon_gs, labels, output_dir)

    return shannon_list, binary_metrics_list, l1norm_list, weighted_unifrac_list


def reformat_pandas(labels, shannon_list, binary_metrics_list, l1norm_list, weighted_unifrac_list):
    """Reformats metrics data into one unified pandas DataFrame.
    
    Parameters
    ----------
    labels : [str]
        List of str for tool names.
    binary_metrics_list : ?
    l1norm_list : ?
    weighted_unifrac_list : ?
    
    Returns
    -------
    Pandas.DataFrame with following columns: metric, rank, tool, value
    """
    # convert unifrac
    pd_weighted_unifrac_list = pd.DataFrame(index=labels, data=weighted_unifrac_list, columns=[c.UNIFRAC, c.UNW_UNIFRAC])
    pd_weighted_unifrac_list = pd_weighted_unifrac_list.stack().reset_index()
    pd_weighted_unifrac_list.columns = ['tool', 'metric', 'value']
    pd_weighted_unifrac_list['rank'] = np.nan

    # convert l1norm
    pd_l1norm_list = pd.DataFrame(l1norm_list, index=labels).stack().reset_index()
    pd_l1norm_list.columns = ['tool', 'rank', 'value']
    pd_l1norm_list['metric'] = c.L1NORM

    # convert binary metrics
    pd_binary_metrics_list = []
    for i, tool in enumerate(labels):
        x = pd.DataFrame([binary_metrics_list[i][rank].get_pretty_dict() for rank in binary_metrics_list[i].keys()]).set_index('rank').stack().reset_index().rename(columns={'level_1': 'metric', 0: 'value'})
        x['tool'] = tool
        pd_binary_metrics_list.append(x)
    pd_binary_metrics_list = pd.concat(pd_binary_metrics_list)
    pd_binary_metrics_list['metric'].replace(['fp',
                                              'tp',
                                              "fn",
                                              "jaccard",
                                              "precision",
                                              "recall"],
                                             [c.FP,
                                              c.TP,
                                              c.FN,
                                              c.JACCARD,
                                              c.PRECISION,
                                              c.RECALL], inplace=True)

    # convert Shannon
    pd_shannon_list = []
    for i, tool in enumerate(labels):
        x = pd.DataFrame([shannon_list[i][rank].get_pretty_dict() for rank in shannon_list[i].keys()]).set_index('rank').stack().reset_index().rename(columns={'level_1': 'metric', 0: 'value'})
        x['tool'] = tool
        pd_shannon_list.append(x)
    pd_shannon_list = pd.concat(pd_shannon_list)
    pd_shannon_list['metric'].replace(['diversity', 'equitability'], [c.SHANNON_DIVERSITY, c.SHANNON_EQUIT], inplace=True)

    # combine metric
    pd_metrics = pd.concat([pd_weighted_unifrac_list, pd_l1norm_list, pd_binary_metrics_list, pd_shannon_list])

    return pd_metrics


def highscore_table(metrics, useranks=['phylum', 'class', 'order', 'family', 'genus']):
    """Compile a ranking table like Figure 3c of CAMI publication.
    
    Note that Figure 3c took into account mean scores for all samples of one of the three
    complexity levels, i.e. 1 for low, 2 for medium, 5 for high.
    Here, I assume that we might be able to do that later, but for now I set "complexity"
    to "dummy".
    
    Parameters
    ----------
    metrics : pd.DataFrame
        Information about metrics of tool performance.
        Must contain columns: metric, rank, tool, value
    useranks : [str]
        Default: 'phylum', 'class', 'order', 'family', 'genus'
        Which ranks should be considered for rank dependent metrics.
        Here we decided to exclude e.g. species, because most profilers
        fail at that rank and we don't want to emphasize on this rank.
    Returns
    -------
    Pandas.DataFrame holding a high scoring table as in Figure 3c.
    """
    # assume for a second that we could handle multiple goldstandard profiles
    pd_metrics = metrics.copy()
    pd_metrics['complexity'] = 'dummy'
    pd_metrics.loc[pd_metrics[pd.isnull(pd_metrics['rank'])].index, 'rank'] = 'rank independent'

    sort_ascendingly = {'L1 norm error': True, 'Unweighted Unifrac error': True,
                        'Recall': False, 'Precision': False}

    # collecting rank scores
    posresults = []
    for (metric, complexity, rank), g in pd_metrics.groupby(['metric', 'complexity', 'rank']):
        if metric in sort_ascendingly:
            if ((rank in useranks) and (metric != 'Unweighted Unifrac error')) or ((rank == 'rank independent') and (metric == 'Unweighted Unifrac error')):
                res = g.groupby('tool').sum().sort_values('value', ascending=sort_ascendingly[metric])
                worstpos = res.shape[0]+1
                res['position'] = range(0, worstpos-1)
                for m in set(pd_metrics['tool'].unique()) - set(res.index):
                    res.loc[m, 'position'] = worstpos
                res['metric'] = metric
                res['complexity'] = complexity
                res['rank'] = rank
                posresults.append(res)
    posresults = pd.concat(posresults)

    # reformat like Figure 3c
    os = []
    for metric, g in posresults.groupby('metric'):
        highscore = g.groupby('tool')['position'].sum().sort_values()
        os.append(pd.DataFrame(["%s (%i)" % (idx, pos) for idx, pos in highscore.iteritems()], columns=[metric]))

    # return reformatted table
    return pd.concat(os, axis=1)  #.T.loc[['Recall', 'Precision', 'L1 norm error', 'Unweighted Unifrac error'],:]


def main():
    parser = argparse.ArgumentParser(description="Compute all metrics for one or more taxonomic profiles")
    parser.add_argument("-g", "--gold_standard_file", help="Gold standard file", required=True)
    parser.add_argument("profiles_files", nargs='+', help="Files of profiles")
    parser.add_argument('-l', '--labels', help="Comma-separated profiles names", required=False)
    parser.add_argument('-o', '--output_dir', help="Directory to write the results to", required=True)
    parser.add_argument('-r', '--by_rank', help="Create a results file per rank", action='store_true')
    args = parser.parse_args()
    labels = get_labels(args.labels, args.profiles_files)
    output_dir = os.path.abspath(args.output_dir)
    make_sure_path_exists(output_dir)
    shannon_list, binary_metrics_list, l1norm_list, weighted_unifrac_list = evaluate(args.gold_standard_file,
                                                                                     args.profiles_files,
                                                                                     labels,
                                                                                     output_dir)
    pd_metrics = reformat_pandas(labels, shannon_list, binary_metrics_list, l1norm_list, weighted_unifrac_list)
    print_by_tool(output_dir, pd_metrics)
    if args.by_rank:
        print_by_rank(output_dir, labels, pd_metrics)


if __name__ == "__main__":
    main()
