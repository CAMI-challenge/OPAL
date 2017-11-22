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
import braycurtis as bc
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
    make_sure_path_exists(os.path.join(output_dir, "by_rank"))
    for rank in c.ALL_RANKS:
        # define ordering of rows, which is given my order of tool labels
        order_rows = labels
        # define ordering of columns, hard coded
        order_columns = [c.UNIFRAC, c.UNW_UNIFRAC, c.L1NORM, c.RECALL, c.PRECISION, c.F1_SCORE, c.TP, c.FP, c.FN, c.JACCARD, c.SHANNON_DIVERSITY, c.SHANNON_EQUIT]
        # subset to those information that either belong to the given rank or are rank independent, i.e. are unifrac values
        table = pd_metrics[(pd_metrics['rank'] == rank) | (pd_metrics['metric'].isin([c.UNIFRAC, c.UNW_UNIFRAC]))]
        # reformat the table with a pivot_table
        table = table.pivot_table(index='tool', columns='metric', values='value')
        # order table
        table = table.loc[order_rows, order_columns]
        # replace np.NaN with string "na" and write resulting table into a file
        table.fillna('na').to_csv(os.path.join(output_dir, "by_rank", rank + ".tsv"), sep='\t')


def print_by_tool(output_dir, pd_metrics):
    make_sure_path_exists(os.path.join(output_dir, "by_tool"))
    # define ordering of columns, hard coded
    order_columns = [c.UNIFRAC, c.UNW_UNIFRAC, c.L1NORM, c.RECALL, c.PRECISION, c.F1_SCORE, c.TP, c.FP, c.FN, c.JACCARD, c.SHANNON_DIVERSITY, c.SHANNON_EQUIT]
    for toolname, pd_metrics_tool in pd_metrics.groupby('tool'):
        table = pd_metrics_tool.pivot_table(index='rank', columns='metric', values='value')
        # little hack to carry unifrac over to every rank
        for unifrac_col in order_columns[:2]:
            table[unifrac_col] = pd_metrics_tool[pd_metrics_tool['metric'] == unifrac_col]['value'].values[0]
        # order table
        table = table.loc[c.ALL_RANKS, order_columns]
        # replace np.NaN with string "na" and write resulting table into a file
        table.fillna('na').to_csv(os.path.join(output_dir, "by_tool", toolname + ".tsv"), sep='\t')


def evaluate(gold_standard_file, profiles_files, labels):
    sample_id, sample_metadata, profile = load_data.open_profile(gold_standard_file)[0]
    gs_rank_to_taxid_to_percentage = load_data.get_rank_to_taxid_to_percentage(profile)
    gs_profile = PF.Profile(sample_metadata=sample_metadata, profile=profile)
    pd_metrics = pd.DataFrame()

    for profile_file, label in zip(profiles_files, labels):
        samples_list = load_data.open_profile(profile_file)
        for sample in samples_list:
            sample_id, sample_metadata, profile = sample

            rank_to_taxid_to_percentage = load_data.get_rank_to_taxid_to_percentage(profile)

            # Unifrac
            pf_profile = PF.Profile(sample_metadata=sample_metadata, profile=profile)
            unifrac = uf.compute_unifrac(gs_profile, pf_profile)

            # Shannon
            shannon = sh.compute_shannon_index(rank_to_taxid_to_percentage)

            # L1 Norm
            l1norm = l1.compute_l1norm(gs_rank_to_taxid_to_percentage, rank_to_taxid_to_percentage)

            # Binary metrics
            binary_metrics = compute_binary_metrics(rank_to_taxid_to_percentage, gs_rank_to_taxid_to_percentage)

            # Bray-Curtis
            braycurtis = bc.braycurtis(gs_rank_to_taxid_to_percentage, rank_to_taxid_to_percentage)

            pd_metrics = pd.concat([pd_metrics, reformat_pandas(sample_id, label, braycurtis, shannon, binary_metrics, l1norm, unifrac)], ignore_index=True)

    return pd_metrics


def reformat_pandas(sample_id, label, braycurtis, shannon, binary_metrics, l1norm, unifrac):
    """Reformats metrics data into one unified pandas DataFrame.
    
    Parameters
    ----------
    sample_id : str
    label : str
        str for tool name.
    braycurtis : float
    shannon : ?
    binary_metrics : ?
    l1norm : float
    unifrac : ?
    
    Returns
    -------
    Pandas.DataFrame with following columns: metric, rank, tool, value
    """
    # convert Unifrac
    pd_unifrac = pd.DataFrame(index=[sample_id], data=[unifrac], columns=[c.UNIFRAC, c.UNW_UNIFRAC]).stack().reset_index()
    pd_unifrac.columns = ['sample', 'metric', 'value']
    pd_unifrac['rank'] = np.nan
    pd_unifrac['tool'] = label

    # convert Shannon
    pd_shannon = pd.DataFrame([shannon[rank].get_pretty_dict() for rank in shannon.keys()]).set_index('rank').stack().reset_index().rename(columns={'level_1': 'metric', 0: 'value'})
    pd_shannon['metric'].replace(['diversity', 'equitability'], [c.SHANNON_DIVERSITY, c.SHANNON_EQUIT], inplace=True)
    pd_shannon['sample'] = sample_id
    pd_shannon['tool'] = label

    # convert L1 Norm
    pd_l1norm = pd.DataFrame(index=[sample_id], data=l1norm).stack().reset_index()
    pd_l1norm.columns = ['sample', 'rank', 'value']
    pd_l1norm['tool'] = label
    pd_l1norm['metric'] = c.L1NORM

    # convert Binary metrics
    pd_binary_metrics = pd.DataFrame([binary_metrics[rank].get_pretty_dict() for rank in binary_metrics.keys()]).set_index('rank').stack().reset_index().rename(columns={'level_1': 'metric', 0: 'value'})
    pd_binary_metrics['metric'].replace(['fp',
                                         'tp',
                                         "fn",
                                         "jaccard",
                                         "precision",
                                         "recall",
                                         "f1"],
                                        [c.FP,
                                         c.TP,
                                         c.FN,
                                         c.JACCARD,
                                         c.PRECISION,
                                         c.RECALL,
                                         c.F1_SCORE], inplace=True)
    pd_binary_metrics['sample'] = sample_id
    pd_binary_metrics['tool'] = label

    # convert Bray-Curtis
    pd_braycurtis = pd.DataFrame(index=[sample_id], data=braycurtis).stack().reset_index()
    pd_braycurtis.columns = ['sample', 'rank', 'value']
    pd_braycurtis['tool'] = label
    pd_braycurtis['metric'] = c.BRAY_CURTIS

    return pd.concat([pd_braycurtis, pd_shannon, pd_binary_metrics, pd_l1norm, pd_unifrac], ignore_index=True)


def plot(pd_metrics, labels, output_dir):
    metrics = [c.UNIFRAC, c.L1NORM, c.RECALL, c.PRECISION, c.FP]
    rank_to_metric_to_toolvalues = defaultdict(dict)

    for rank in c.PHYLUM_SPECIES:
        for metric in metrics:
            rank_to_metric_to_toolvalues[rank][metric] = []
        table1 = pd_metrics[(pd_metrics['rank'] == rank) | (pd_metrics['metric'].isin([c.UNIFRAC]))]
        max_fp = table1[table1['metric'] == c.FP]['value'].max()

        for label in labels:
            table2 = table1[table1['tool'] == label]

            unifrac = table2[table2['metric'] == c.UNIFRAC]['value'].values[0]
            rank_to_metric_to_toolvalues[rank][c.UNIFRAC].append(unifrac / 16)

            l1norm = table2[table2['metric'] == c.L1NORM]['value'].values
            rank_to_metric_to_toolvalues[rank][c.L1NORM].append(l1norm[0] / 2.0 if len(l1norm) > 0 else None)

            recall = table2[table2['metric'] == c.RECALL]['value'].values
            rank_to_metric_to_toolvalues[rank][c.RECALL].append(recall[0] if len(recall) > 0 else None)

            precision = table2[table2['metric'] == c.PRECISION]['value'].values
            rank_to_metric_to_toolvalues[rank][c.PRECISION].append(precision[0] if len(precision) > 0 else None)

            fp = table2[table2['metric'] == c.FP]['value'].values
            if max_fp > 0:
                rank_to_metric_to_toolvalues[rank][c.FP].append(fp[0] / max_fp if len(fp) > 0 else None)
            else:
                rank_to_metric_to_toolvalues[rank][c.FP].append(fp[0] if len(fp) > 0 else None)

    pl.spider_plot(metrics,
                   labels,
                   rank_to_metric_to_toolvalues,
                   output_dir,
                   'spider_plot',
                   ['b', 'g', 'r', 'k', 'm'])

    pl.spider_plot([c.RECALL, c.PRECISION],
                   labels,
                   rank_to_metric_to_toolvalues,
                   output_dir,
                   'spider_plot_recall_precision',
                   ['r', 'k'],
                   grid_points=[0.2, 0.4, 0.6, 0.8, 1.0],
                   fill=True)

    # TODO: re-enable
    # rank_to_shannon_gs = sh.compute_shannon_index(gs_rank_to_taxid_to_percentage)
    # pl.plot_shannon(shannon_list, rank_to_shannon_gs, labels, output_dir)

    ## pl.plot_braycurtis_l1norm(braycurtis_list, l1norm_list, labels, output_dir)


def main():
    parser = argparse.ArgumentParser(description="Compute all metrics for one or more taxonomic profiles")
    parser.add_argument("-g", "--gold_standard_file", help="Gold standard file", required=True)
    parser.add_argument("profiles_files", nargs='+', help="Files of profiles")
    parser.add_argument('-l', '--labels', help="Comma-separated profiles names", required=False)
    parser.add_argument('-o', '--output_dir', help="Directory to write the results to", required=True)
    args = parser.parse_args()
    labels = get_labels(args.labels, args.profiles_files)
    output_dir = os.path.abspath(args.output_dir)
    make_sure_path_exists(output_dir)
    pd_metrics = evaluate(args.gold_standard_file,
                          args.profiles_files,
                          labels)
    pd_metrics[['tool', 'rank', 'metric', 'sample', 'value']].fillna('na').to_csv(os.path.join(output_dir, "results.tsv"), sep='\t', index=False)
    print_by_tool(output_dir, pd_metrics)
    print_by_rank(output_dir, labels, pd_metrics)
    plot(pd_metrics, labels, output_dir)


if __name__ == "__main__":
    main()
