#!/usr/bin/env python3

import os
import sys
import errno
import argparse
import os.path
import pandas as pd
import numpy as np
import logging
import shlex
from cami_opal import evaluate
from cami_opal import rankings as rk
from cami_opal import html_opal as html
from cami_opal import plots as pl
from cami_opal.utils import load_data
from cami_opal.utils import constants as c
from version import __version__


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
            logging.getLogger('opal').critical('The number of labels does not match the number of files of profiles. Please check parameter -l, --labels.')
            exit(1)
        return labels_list
    tool_id = []
    for profile_file in profiles_files:
        tool_id.append(profile_file.split('/')[-1])
    return tool_id


def get_time_memory(time, memory, profiles_files):
    time_list = []
    memory_list = []
    if time:
        time_list = [float(x.strip()) for x in time.split(',')]
        if len(time_list) != len(profiles_files):
            logging.getLogger('opal').critical('The number of running times does not match the number of files of profiles. Please check parameter --time.')
            exit(1)
    if memory:
        memory_list = [float(x.strip()) for x in memory.split(',')]
        if len(memory_list) != len(profiles_files):
            logging.getLogger('opal').critical('The number of memory usages does not match the number of files of profiles. Please check parameter --memory.')
            exit(1)
    return time_list, memory_list


def print_by_rank(output_dir, labels, pd_metrics):
    make_sure_path_exists(os.path.join(output_dir, "by_rank"))
    # define ordering of rows, which is given my order of tool labels
    order_rows = labels
    # define ordering of columns, hard coded
    order_columns = [c.SUM_ABUNDANCES, c.UNIFRAC, c.UNW_UNIFRAC, c.L1NORM, c.RECALL, c.PRECISION, c.F1_SCORE, c.TP, c.FP, c.FN, c.OTUS, c.JACCARD, c.SHANNON_DIVERSITY, c.SHANNON_EQUIT, c.BRAY_CURTIS]
    if c.FP + c.UNFILTERED_SUF in pd_metrics['metric'].values:
        order_columns += [metric + c.UNFILTERED_SUF for metric in order_columns]
    for rank in c.ALL_RANKS:
        order_columns_rank = order_columns
        # subset to information that either belongs to the given rank or is rank independent, i.e. are unifrac values
        table = pd_metrics[(pd_metrics['rank'] == rank) | (pd_metrics['metric'].isin([c.UNIFRAC, c.UNW_UNIFRAC, c.UNIFRAC + c.UNFILTERED_SUF, c.UNW_UNIFRAC + c.UNFILTERED_SUF]))]
        # reformat the table with a pivot_table
        table = table.pivot_table(index=['tool', 'sample'], columns='metric', values='value')
        if len(table.columns) < len(order_columns):
            order_columns_rank = [x for x in order_columns if x in table.columns]
        # select only tools in labels and get rid of gold standard
        table = table.loc[pd.IndexSlice[order_rows, :], order_columns_rank]
        # define categorical column for ordering rows by tools
        table['tool_cat'] = pd.Categorical(table.index.get_level_values('tool'), categories=order_rows, ordered=True)
        # order table
        table = table.sort_values('tool_cat')
        table = table.loc[:, order_columns_rank]
        # replace np.NaN with string "na" and write resulting table into a file
        table.fillna('na').to_csv(os.path.join(output_dir, "by_rank", rank + ".tsv"), sep='\t')


def print_by_tool(output_dir, pd_metrics):
    make_sure_path_exists(os.path.join(output_dir, "by_tool"))
    # define ordering of columns, hard coded
    order_columns = [c.UNIFRAC, c.UNW_UNIFRAC, c.L1NORM, c.RECALL, c.PRECISION, c.F1_SCORE, c.TP, c.FP, c.FN, c.OTUS, c.JACCARD, c.SHANNON_DIVERSITY, c.SHANNON_EQUIT, c.BRAY_CURTIS]
    unifrac_list = [c.UNIFRAC, c.UNW_UNIFRAC]
    if c.FP + c.UNFILTERED_SUF in pd_metrics['metric'].values:
        order_columns += [metric + c.UNFILTERED_SUF for metric in order_columns]
        unifrac_list += [c.UNIFRAC + c.UNFILTERED_SUF, c.UNW_UNIFRAC + c.UNFILTERED_SUF]
    for toolname, pd_metrics_tool in pd_metrics.groupby('tool'):
        if toolname == c.GS:
            continue
        table = pd_metrics_tool.pivot_table(index=['rank', 'sample'], columns='metric', values='value')
        # little hack to carry unifrac over to every rank
        for unifrac_col in unifrac_list:
            table[unifrac_col] = pd_metrics_tool[pd_metrics_tool['metric'] == unifrac_col]['value'].values[0]
        # order table
        table['rank_cat'] = pd.Categorical(table.index.get_level_values('rank'), categories=c.ALL_RANKS, ordered=True)
        table = table.sort_values('rank_cat')
        table = table.loc[:, order_columns]
        # replace np.NaN with string "na" and write resulting table into a file
        table.fillna('na').to_csv(os.path.join(output_dir, "by_tool", toolname + ".tsv"), sep='\t')


def create_output_directories(output_dir, labels):
    make_sure_path_exists(os.path.join(output_dir, 'gold_standard'))
    for label in labels:
        make_sure_path_exists(os.path.join(output_dir, "by_tool", label.replace(' ', '_')))


def concat_pd(labels, metric, values, pd_metrics):
    df = pd.DataFrame({'tool': labels, 'value': values})
    df['sample'] = np.nan
    df['metric'] = metric
    df['rank'] = np.nan
    return pd.concat([pd_metrics, df], ignore_index=True, sort=False)


def concat_time_memory(labels, time_list, memory_list, pd_metrics):
    if time_list:
        pd_metrics = concat_pd(labels, 'time', time_list, pd_metrics)
    if memory_list:
        pd_metrics = concat_pd(labels, 'memory', memory_list, pd_metrics)
    return pd_metrics


def get_logger(output_dir, silent):
    logger = logging.getLogger('opal')
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    logging_fh = logging.FileHandler(os.path.join(output_dir, 'log.txt'))
    logging_fh.setFormatter(formatter)
    logger.addHandler(logging_fh)

    logger.info(' '.join(map(shlex.quote, sys.argv)))

    if not silent:
        logging_stdout = logging.StreamHandler(sys.stdout)
        logging_stdout.setFormatter(formatter)
        logger.addHandler(logging_stdout)
    return logger


def main():
    parser = argparse.ArgumentParser(description='OPAL: Open-community Profiling Assessment tooL', add_help=False)
    group1 = parser.add_argument_group('required arguments')
    group1.add_argument('profiles_files', nargs='+', help='Files of profiles')
    group1.add_argument('-g', '--gold_standard_file', help='Gold standard file', required=True)
    group1.add_argument('-o', '--output_dir', help='Directory to write the results to', required=True)

    group2 = parser.add_argument_group('optional arguments')
    group2.add_argument('-n', '--normalize', help='Normalize samples', action='store_true')
    group2.add_argument('-f', '--filter', help='Filter out the predictions with the smallest relative abundances summing up to [FILTER]%% within a rank', type=float)
    group2.add_argument('-p', '--plot_abundances', help='Plot abundances in the gold standard (can take some minutes)', action='store_true')
    group2.add_argument('-l', '--labels', help='Comma-separated profiles names', required=False)
    group2.add_argument('-t', '--time', help='Comma-separated runtimes in hours', required=False)
    group2.add_argument('-m', '--memory', help='Comma-separated memory usages in gigabytes', required=False)
    group2.add_argument('-d', '--desc', help='Description for HTML page', required=False)
    group2.add_argument('-r', '--ranks', help='Highest and lowest taxonomic ranks to consider in performance rankings, comma-separated. Valid ranks: superkingdom, phylum, class, order, family, genus, species, strain (default:superkingdom,species)', required=False)
    group2.add_argument('--metrics_plot_rel', help='Metrics for spider plot of relative performances, first character, comma-separated. Valid metrics: w:weighted Unifrac, l:L1 norm, c:completeness, p:purity, f:false positives, t:true positives (default: w,l,c,p,f)', required=False)
    group2.add_argument('--metrics_plot_abs', help='Metrics for spider plot of absolute performances, first character, comma-separated. Valid metrics: c:completeness, p:purity, b:Bray-Curtis (default: c,p)', required=False)
    group2.add_argument('--silent', help='Silent mode', action='store_true')
    group2.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    group2.add_argument('-h', '--help', action='help', help='Show this help message and exit')

    group3 = parser.add_argument_group('UniFrac arguments')
    group3.add_argument('-b', '--branch_length_function', help='UniFrac tree branch length function (default: "lambda x: 1/x", where x=tree depth)', required=False, default='lambda x: 1/x')
    group3.add_argument('--normalized_unifrac', help='Compute normalized version of weighted UniFrac by dividing by the theoretical max unweighted UniFrac', action='store_true')

    args = parser.parse_args()
    output_dir = os.path.abspath(args.output_dir)
    make_sure_path_exists(output_dir)
    labels = get_labels(args.labels, args.profiles_files)

    create_output_directories(output_dir, labels)

    logger = get_logger(args.output_dir, args.silent)

    logger.info('Loading profiles...')
    sample_ids_list, gs_samples_list, profiles_list_to_samples_list = load_data.load_profiles(args.gold_standard_file,
                                                                                    args.profiles_files,
                                                                                    args.normalize)
    logger.info('done')

    plots_list = []
    if args.plot_abundances:
        logger.info('Plotting gold standard abundances...')
        plots_list += pl.plot_samples_hist(gs_samples_list, sample_ids_list, output_dir)
        logger.info('done')

    logger.info('Computing metrics...')
    pd_metrics = evaluate.evaluate_main(gs_samples_list,
                                        profiles_list_to_samples_list,
                                        labels,
                                        args.filter,
                                        args.branch_length_function,
                                        args.normalized_unifrac)
    time_list, memory_list = get_time_memory(args.time, args.memory, args.profiles_files)
    if time_list or memory_list:
        pd_metrics = concat_time_memory(labels, time_list, memory_list, pd_metrics)
    logger.info('done')

    logger.info('Saving computed metrics...')
    pd_metrics[['tool', 'rank', 'metric', 'sample', 'value']].fillna('na').to_csv(os.path.join(output_dir, 'results.tsv'), sep='\t', index=False)
    print_by_tool(output_dir, pd_metrics)
    print_by_rank(output_dir, labels, pd_metrics)
    logger.info('done')

    logger.info('Creating beta diversity plots...')
    plots_list += pl.plot_beta_diversity(gs_samples_list, profiles_list_to_samples_list, sample_ids_list, labels, output_dir)
    logger.info('done')

    logger.info('Creating rarefaction curves...')
    plots_list += pl.plot_rarefaction_curves(gs_samples_list, output_dir)
    plots_list += pl.plot_rarefaction_curves(gs_samples_list, output_dir, log_scale=True)
    logger.info('done')

    logger.info('Creating more plots...')
    plots_list += pl.plot_all(pd_metrics, labels, output_dir, args.metrics_plot_rel, args.metrics_plot_abs)
    logger.info('done')

    logger.info('Computing rankings...')
    pd_rankings, ranks_scored = rk.highscore_table(pd_metrics, args.ranks)
    logger.info('done')

    if time_list or memory_list:
        logger.info('Plotting computing efficiency...')
        plots_list += pl.plot_time_memory(time_list, memory_list, labels, output_dir)
        logger.info('done')

    logger.info('Creating HTML page...')
    html.create_html(pd_rankings, ranks_scored, pd_metrics, labels, sample_ids_list, plots_list, output_dir, args.desc)
    logger.info('done')

    logger.info('OPAL finished successfully. All results have been saved to {}'.format(output_dir))


if __name__ == "__main__":
    main()
