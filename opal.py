#!/usr/bin/env python3

import os
import sys
import errno
import argparse
import os.path
import pandas as pd
import numpy as np
import logging
from src import l1norm as l1
from src import binary_metrics as bm
from src import unifrac_distance as uf
from src import shannon as sh
from src import braycurtis as bc
from src import rankings as rk
from src import html_opal as html
from src import plots as pl
from src.utils import load_data
from src.utils import ProfilingTools as PF
from src.utils import constants as c
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


def compute_binary_metrics(query_profile, query_truth):
    all_metrics = bm.compute_tree_metrics(query_profile, query_truth)
    return all_metrics


def print_by_rank(output_dir, labels, pd_metrics):
    make_sure_path_exists(os.path.join(output_dir, "by_rank"))
    for rank in c.ALL_RANKS:
        # define ordering of rows, which is given my order of tool labels
        order_rows = labels
        # define ordering of columns, hard coded
        order_columns = [c.UNIFRAC, c.UNW_UNIFRAC, c.L1NORM, c.RECALL, c.PRECISION, c.F1_SCORE, c.TP, c.FP, c.FN, c.OTUS, c.JACCARD, c.SHANNON_DIVERSITY, c.SHANNON_EQUIT, c.BRAY_CURTIS]
        # subset to those information that either belong to the given rank or are rank independent, i.e. are unifrac values
        table = pd_metrics[(pd_metrics['rank'] == rank) | (pd_metrics['metric'].isin([c.UNIFRAC, c.UNW_UNIFRAC]))]
        # reformat the table with a pivot_table
        table = table.pivot_table(index=['tool', 'sample'], columns='metric', values='value')
        # select only tools in labels and get rid of gold standard
        table = table.loc[pd.IndexSlice[order_rows,:], order_columns]
        # define categorical column for ordering rows by tools
        table['tool_cat'] = pd.Categorical(table.index.get_level_values('tool'), categories=order_rows, ordered=True)
        # order table
        table = table.sort_values('tool_cat')
        table = table.loc[:, order_columns]
        # replace np.NaN with string "na" and write resulting table into a file
        table.fillna('na').to_csv(os.path.join(output_dir, "by_rank", rank + ".tsv"), sep='\t')


def print_by_tool(output_dir, pd_metrics):
    make_sure_path_exists(os.path.join(output_dir, "by_tool"))
    # define ordering of columns, hard coded
    order_columns = [c.UNIFRAC, c.UNW_UNIFRAC, c.L1NORM, c.RECALL, c.PRECISION, c.F1_SCORE, c.TP, c.FP, c.FN, c.OTUS, c.JACCARD, c.SHANNON_DIVERSITY, c.SHANNON_EQUIT, c.BRAY_CURTIS]
    for toolname, pd_metrics_tool in pd_metrics.groupby('tool'):
        if toolname == c.GS:
            continue
        table = pd_metrics_tool.pivot_table(index=['rank', 'sample'], columns='metric', values='value')
        # little hack to carry unifrac over to every rank
        for unifrac_col in order_columns[:2]:
            table[unifrac_col] = pd_metrics_tool[pd_metrics_tool['metric'] == unifrac_col]['value'].values[0]
        # order table
        table['rank_cat'] = pd.Categorical(table.index.get_level_values('rank'), categories=c.ALL_RANKS, ordered=True)
        table = table.sort_values('rank_cat')
        table = table.loc[:, order_columns]
        # replace np.NaN with string "na" and write resulting table into a file
        table.fillna('na').to_csv(os.path.join(output_dir, "by_tool", toolname + ".tsv"), sep='\t')


def compute_metrics(sample_metadata, profile, gs_pf_profile, gs_rank_to_taxid_to_percentage, rank_to_taxid_to_percentage):
    # Unifrac
    if isinstance(profile, PF.Profile):
        pf_profile = profile
    else:
        pf_profile = PF.Profile(sample_metadata=sample_metadata, profile=profile)
    unifrac = uf.compute_unifrac(gs_pf_profile, pf_profile)

    # Shannon
    shannon = sh.compute_shannon_index(rank_to_taxid_to_percentage)

    # L1 Norm
    l1norm = l1.compute_l1norm(gs_rank_to_taxid_to_percentage, rank_to_taxid_to_percentage)

    # Binary metrics
    binary_metrics = compute_binary_metrics(rank_to_taxid_to_percentage, gs_rank_to_taxid_to_percentage)

    # Bray-Curtis
    braycurtis = bc.braycurtis(gs_rank_to_taxid_to_percentage, rank_to_taxid_to_percentage)

    return unifrac, shannon, l1norm, binary_metrics, braycurtis


def load_profiles(gold_standard_file, profiles_files, no_normalization):
    normalize = False if no_normalization else True

    gs_samples_list = load_data.open_profile(gold_standard_file, normalize)
    sample_ids_list = []
    for sample in gs_samples_list:
        sample_id, sample_metadata, profile = sample
        sample_ids_list.append(sample_id)

    profiles_list_to_samples_list = []
    for profile_file in profiles_files:
        profiles_list_to_samples_list.append(load_data.open_profile(profile_file, normalize))

    return sample_ids_list, gs_samples_list, profiles_list_to_samples_list


def evaluate(gs_samples_list, profiles_list_to_samples_list, labels):
    gs_id_to_rank_to_taxid_to_percentage = {}
    gs_id_to_pf_profile = {}
    pd_metrics = pd.DataFrame()

    for sample in gs_samples_list:
        sample_id, sample_metadata, profile = sample
        gs_id_to_rank_to_taxid_to_percentage[sample_id] = load_data.get_rank_to_taxid_to_percentage(profile)
        gs_id_to_pf_profile[sample_id] = PF.Profile(sample_metadata=sample_metadata, profile=profile)
        unifrac, shannon, l1norm, binary_metrics, braycurtis = compute_metrics(sample_metadata,
                                                                               gs_id_to_pf_profile[sample_id],
                                                                               gs_id_to_pf_profile[sample_id],
                                                                               gs_id_to_rank_to_taxid_to_percentage[sample_id],
                                                                               gs_id_to_rank_to_taxid_to_percentage[sample_id])
        pd_metrics = pd.concat([pd_metrics, reformat_pandas(sample_id, c.GS, braycurtis, shannon, binary_metrics, l1norm, unifrac)], ignore_index=True)

    one_profile_assessed = False
    for samples_list, label in zip(profiles_list_to_samples_list, labels):
        for sample in samples_list:
            sample_id, sample_metadata, profile = sample

            # match the sample id of the gold standard and the predictions
            if sample_id in gs_id_to_rank_to_taxid_to_percentage:
                gs_rank_to_taxid_to_percentage = gs_id_to_rank_to_taxid_to_percentage[sample_id]
                gs_pf_profile = gs_id_to_pf_profile[sample_id]
            else:
                sys.stderr.write("Skipping assessment of {} for sample {}. Make sure the SampleID of the gold standard and the profile are identical.\n".format(label, sample_id))
                continue

            rank_to_taxid_to_percentage = load_data.get_rank_to_taxid_to_percentage(profile)

            unifrac, shannon, l1norm, binary_metrics, braycurtis = compute_metrics(sample_metadata, profile, gs_pf_profile,
                                                                                   gs_rank_to_taxid_to_percentage,
                                                                                   rank_to_taxid_to_percentage)
            pd_metrics = pd.concat([pd_metrics, reformat_pandas(sample_id, label, braycurtis, shannon, binary_metrics, l1norm, unifrac)], ignore_index=True)
            one_profile_assessed = True

    if not one_profile_assessed:
        logging.getLogger('opal').critical("No profile could be evaluated.")
        exit(1)

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
                                         "f1",
                                         "otus"],
                                        [c.FP,
                                         c.TP,
                                         c.FN,
                                         c.JACCARD,
                                         c.PRECISION,
                                         c.RECALL,
                                         c.F1_SCORE,
                                         c.OTUS], inplace=True)
    pd_binary_metrics['sample'] = sample_id
    pd_binary_metrics['tool'] = label

    # convert Bray-Curtis
    pd_braycurtis = pd.DataFrame(index=[sample_id], data=braycurtis).stack().reset_index()
    pd_braycurtis.columns = ['sample', 'rank', 'value']
    pd_braycurtis['tool'] = label
    pd_braycurtis['metric'] = c.BRAY_CURTIS

    return pd.concat([pd_braycurtis, pd_shannon, pd_binary_metrics, pd_l1norm, pd_unifrac], ignore_index=True, sort=False)


def create_output_directories(output_dir, labels):
    make_sure_path_exists(os.path.join(output_dir, 'gold_standard'))
    for label in labels:
        make_sure_path_exists(os.path.join(output_dir, "by_tool", label))


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
    group2.add_argument('-n', '--no_normalization', help='Do not normalize samples', action='store_true')
    group2.add_argument('-p', '--plot_abundances', help='Plot abundances in the gold standard (can take some minutes)', action='store_true')
    group2.add_argument('-l', '--labels', help='Comma-separated profiles names', required=False)
    group2.add_argument('-t', '--time', help='Comma-separated runtimes in hours', required=False)
    group2.add_argument('-m', '--memory', help='Comma-separated memory usages in gigabytes', required=False)
    group2.add_argument('-d', '--desc', help='Description for HTML page', required=False)
    group2.add_argument('--silent', help='Silent mode', action='store_true')
    group2.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    group2.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    args = parser.parse_args()
    output_dir = os.path.abspath(args.output_dir)
    make_sure_path_exists(output_dir)
    labels = get_labels(args.labels, args.profiles_files)

    create_output_directories(output_dir, labels)

    logger = get_logger(args.output_dir, args.silent)

    logger.info('Loading profiles...')
    sample_ids_list, gs_samples_list, profiles_list_to_samples_list = load_profiles(args.gold_standard_file,
                                                                                    args.profiles_files,
                                                                                    args.no_normalization)
    logger.info('done')

    plots_list = []
    if args.plot_abundances:
        logger.info('Plotting gold standard abundances...')
        plots_list += pl.plot_samples_hist(gs_samples_list, sample_ids_list, output_dir)
        logger.info('done')

    logger.info('Computing metrics...')
    pd_metrics = evaluate(gs_samples_list,
                          profiles_list_to_samples_list,
                          labels)
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

    logger.info('Creating spider plots...')
    plots_list += pl.plot_all(pd_metrics, labels, output_dir)
    logger.info('done')

    logger.info('Computing rankings...')
    pd_rankings = rk.highscore_table(pd_metrics)
    logger.info('done')

    if time_list or memory_list:
        logger.info('Plotting computing efficiency...')
        plots_list += pl.plot_time_memory(time_list, memory_list, labels, output_dir)
        logger.info('done')

    logger.info('Creating HTML page...')
    html.create_html(pd_rankings, pd_metrics, labels, sample_ids_list, plots_list, output_dir, args.desc)
    logger.info('done')

    logger.info('OPAL finished successfully. All results have been saved to {}'.format(output_dir))


if __name__ == "__main__":
    main()
