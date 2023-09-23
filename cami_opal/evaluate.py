from cami_opal import l1norm as l1
from cami_opal import binary_metrics as bm
from cami_opal import unifrac_distance as uf
from cami_opal import shannon as sh
from cami_opal import braycurtis as bc
from cami_opal.utils import ProfilingTools as PF
from cami_opal.utils import ProfilingToolsCAMI as PFCAMI
from cami_opal.utils import constants as c
from cami_opal.utils import load_data

import pandas as pd
import numpy as np
import logging


def compute_binary_metrics(query_profile, query_truth):
    all_metrics = bm.compute_tree_metrics(query_profile, query_truth)
    return all_metrics


def compute_metrics(sample_metadata, profile, gs_pf_profile, profile_cami, gs_pf_profile_cami,
                    gs_rank_to_taxid_to_percentage, rank_to_taxid_to_percentage, branch_length_fun, normalized_unifrac):
    # Unifrac
    if isinstance(profile, PF.Profile):
        pf_profile = profile
        pf_profile_cami = profile_cami
    else:
        pf_profile = PF.Profile(sample_metadata=sample_metadata, profile=profile, branch_length_fun=branch_length_fun)
        pf_profile_cami = PFCAMI.ProfileCAMI(sample_metadata=sample_metadata, profile=profile_cami)
    unifrac = uf.compute_unifrac(gs_pf_profile, pf_profile, normalized_unifrac)
    unifrac_cami = uf.compute_unifrac(gs_pf_profile_cami, pf_profile_cami)

    # Shannon
    shannon = sh.compute_shannon_index(rank_to_taxid_to_percentage)

    # L1 Norm
    l1norm = l1.compute_l1norm(gs_rank_to_taxid_to_percentage, rank_to_taxid_to_percentage)

    # Binary metrics
    binary_metrics = compute_binary_metrics(rank_to_taxid_to_percentage, gs_rank_to_taxid_to_percentage)

    # Bray-Curtis
    braycurtis = bc.braycurtis(gs_rank_to_taxid_to_percentage, rank_to_taxid_to_percentage)

    # Sum of abundances and taxon counts
    rank_to_sum = {}
    rank_to_ntaxa = {}
    for rank in rank_to_taxid_to_percentage:
        rank_to_sum[rank] = sum(rank_to_taxid_to_percentage[rank].values())
        rank_to_ntaxa[rank] = len(rank_to_taxid_to_percentage[rank])

    return unifrac, unifrac_cami, shannon, l1norm, binary_metrics, braycurtis, rank_to_sum, rank_to_ntaxa


def reformat_pandas(sample_id, label, braycurtis, shannon, binary_metrics, l1norm, unifrac, unifrac_cami, rank_to_sum,
                    rank_to_ntaxa, rename_as_unfiltered=False):
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
    if not l1norm:
        return pd.DataFrame()

    # convert Unifrac
    pd_unifrac = pd.DataFrame(index=[sample_id], data=[unifrac],
                              columns=[c.UNIFRAC, c.UNW_UNIFRAC]).stack().reset_index()
    pd_unifrac.columns = ['sample', 'metric', 'value']
    pd_unifrac['rank'] = np.nan
    pd_unifrac['tool'] = label

    # convert Unifrac CAMI
    pd_unifrac_cami = pd.DataFrame(index=[sample_id], data=[unifrac_cami],
                                   columns=[c.UNIFRAC_CAMI, c.UNW_UNIFRAC_CAMI]).stack().reset_index()
    pd_unifrac_cami.columns = ['sample', 'metric', 'value']
    pd_unifrac_cami['rank'] = np.nan
    pd_unifrac_cami['tool'] = label

    # convert Shannon
    pd_shannon = pd.DataFrame([shannon[rank].get_pretty_dict() for rank in shannon.keys()]).set_index(
        'rank').stack().reset_index().rename(columns={'level_1': 'metric', 0: 'value'})
    pd_shannon['metric'].replace(['diversity', 'equitability'], [c.SHANNON_DIVERSITY, c.SHANNON_EQUIT], inplace=True)
    pd_shannon['sample'] = sample_id
    pd_shannon['tool'] = label

    # convert L1 Norm
    pd_l1norm = pd.DataFrame(index=[sample_id], data=l1norm).stack().reset_index()
    pd_l1norm.columns = ['sample', 'rank', 'value']
    pd_l1norm['tool'] = label
    pd_l1norm['metric'] = c.L1NORM

    # convert Binary metrics
    pd_binary_metrics = pd.DataFrame(
        [binary_metrics[rank].get_pretty_dict() for rank in binary_metrics.keys()]).set_index(
        'rank').stack().reset_index().rename(columns={'level_1': 'metric', 0: 'value'})
    pd_binary_metrics['metric'].replace(['fp', 'tp', 'fn', 'jaccard', 'precision', 'recall', 'f1'],
                                        [c.FP, c.TP, c.FN, c.JACCARD, c.PRECISION, c.RECALL, c.F1_SCORE],
                                        inplace=True)
    pd_binary_metrics['sample'] = sample_id
    pd_binary_metrics['tool'] = label

    # convert Bray-Curtis
    pd_braycurtis = pd.DataFrame(index=[sample_id], data=braycurtis).stack().reset_index()
    pd_braycurtis.columns = ['sample', 'rank', 'value']
    pd_braycurtis['tool'] = label
    pd_braycurtis['metric'] = c.BRAY_CURTIS

    # convert Sum of abundances
    pd_sum = pd.DataFrame(index=[sample_id], data=rank_to_sum).stack().reset_index()
    pd_sum.columns = ['sample', 'rank', 'value']
    pd_sum['tool'] = label
    pd_sum['metric'] = c.SUM_ABUNDANCES
    pd_sum = pd_sum[pd_sum['rank'].isin(c.ALL_RANKS)]
    pd_sum['value'] = pd_sum['value'] / 100.0

    # convert Taxon counts
    pd_ntaxa = pd.DataFrame(index=[sample_id], data=rank_to_ntaxa).stack().reset_index()
    pd_ntaxa.columns = ['sample', 'rank', 'value']
    pd_ntaxa['tool'] = label
    pd_ntaxa['metric'] = c.OTUS
    pd_ntaxa = pd_ntaxa[pd_ntaxa['rank'].isin(c.ALL_RANKS)]
    pd_ntaxa['value'] = pd_ntaxa['value']

    pd_formatted = pd.concat(
        [pd_braycurtis, pd_shannon, pd_binary_metrics, pd_l1norm, pd_unifrac, pd_unifrac_cami, pd_sum, pd_ntaxa],
        ignore_index=True, sort=False)

    if rename_as_unfiltered:
        metrics_list = pd_formatted['metric'].unique().tolist()
        pd_formatted['metric'].replace(metrics_list, [metric + c.UNFILTERED_SUF for metric in metrics_list],
                                       inplace=True)

    return pd_formatted


def get_confusion_df(gs_rank_to_taxid_to_percentage, rank_to_taxid_to_percentage):
    def dict_to_pandas(pct_dict):
        return pd.DataFrame.from_dict(pct_dict).stack().reset_index().dropna(subset=[0]).rename(
            columns={'level_0': 'taxid', 'level_1': 'rank', 0: 'pct'})

    df_gs = dict_to_pandas(gs_rank_to_taxid_to_percentage)
    df_pred = dict_to_pandas(rank_to_taxid_to_percentage)
    df_pred['classification'] = df_pred['taxid'].isin(df_gs['taxid'])
    df_pred['classification'] = np.where(df_pred['classification'], 'TP', 'FP')
    df_pred = pd.merge(df_pred, df_gs, on='taxid', sort=False, how='outer')
    df_pred['classification'] = df_pred['classification'].fillna('FN')
    df_pred['rank_x'] = np.where(df_pred['rank_x'].isna(), df_pred['rank_y'], df_pred['rank_x'])
    df_pred = df_pred.rename(columns={'rank_x': 'rank', 'pct_x': 'pct', 'pct_y': 'pct_gs'}).drop('rank_y', axis=1)
    return df_pred


def evaluate_gs(gs_samples_list, filter_tail_percentage, branch_length_fun, normalized_unifrac,
                gs_id_to_rank_to_taxid_to_percentage, gs_id_to_pf_profile, gs_id_to_pf_profile_cami, skip_gs=False):
    pd_metrics = pd.DataFrame()
    pd_confusion = pd.DataFrame()
    for sample in gs_samples_list:
        sample_id, sample_metadata, profile = sample
        gs_id_to_rank_to_taxid_to_percentage[sample_id] = load_data.get_rank_to_taxid_to_percentage(profile)
        gs_id_to_pf_profile[sample_id] = PF.Profile(sample_metadata=sample_metadata, profile=profile)
        gs_id_to_pf_profile_cami[sample_id] = PFCAMI.ProfileCAMI(sample_metadata=sample_metadata, profile=profile)
        if not skip_gs:
            unifrac, unifrac_cami, shannon, l1norm, binary_metrics, braycurtis, rank_to_sum, rank_to_ntaxa = \
                compute_metrics(sample_metadata,
                                gs_id_to_pf_profile[sample_id], gs_id_to_pf_profile[sample_id],
                                gs_id_to_pf_profile_cami[sample_id], gs_id_to_pf_profile_cami[sample_id],
                                gs_id_to_rank_to_taxid_to_percentage[sample_id],
                                gs_id_to_rank_to_taxid_to_percentage[sample_id],
                                branch_length_fun, normalized_unifrac)
            pd_metrics = pd.concat([pd_metrics,
                                    reformat_pandas(sample_id, c.GS, braycurtis, shannon, binary_metrics, l1norm,
                                                    unifrac, unifrac_cami, rank_to_sum, rank_to_ntaxa)],
                                   ignore_index=True)

            pd_sample_confusion = get_confusion_df(gs_id_to_rank_to_taxid_to_percentage[sample_id], gs_id_to_rank_to_taxid_to_percentage[sample_id])
            pd_sample_confusion['sample'] = sample_id
            pd_sample_confusion['tool'] = c.GS
            pd_confusion = pd.concat([pd_confusion, pd_sample_confusion], ignore_index=True)

    if filter_tail_percentage and not skip_gs:
        metrics_list = pd_metrics['metric'].unique().tolist()
        pd_metrics_copy = pd_metrics.copy()
        pd_metrics_copy['metric'].replace(metrics_list, [metric + c.UNFILTERED_SUF for metric in metrics_list],
                                          inplace=True)
        pd_metrics = pd.concat([pd_metrics, pd_metrics_copy], ignore_index=True)
    return pd_metrics, pd_confusion


def evaluate_main(gs_samples_list, profiles_list_to_samples_list, labels, filter_tail_percentage, branch_length,
                  normalized_unifrac, skip_gs=False):
    gs_id_to_rank_to_taxid_to_percentage = {}
    gs_id_to_pf_profile = {}
    gs_id_to_pf_profile_cami = {}
    branch_length_fun = PF.Profile.get_branch_length_function(branch_length)

    pd_metrics, pd_confusion = evaluate_gs(gs_samples_list, filter_tail_percentage, branch_length_fun, normalized_unifrac,
                             gs_id_to_rank_to_taxid_to_percentage, gs_id_to_pf_profile, gs_id_to_pf_profile_cami,
                             skip_gs)

    one_profile_assessed = False
    for samples_list, label in zip(profiles_list_to_samples_list, labels):
        for sample in samples_list:
            sample_id, sample_metadata, profile = sample

            # match the sample id of the gold standard and the predictions
            if sample_id in gs_id_to_rank_to_taxid_to_percentage:
                gs_rank_to_taxid_to_percentage = gs_id_to_rank_to_taxid_to_percentage[sample_id]
                gs_pf_profile = gs_id_to_pf_profile[sample_id]
                gs_pf_profile_cami = gs_id_to_pf_profile_cami[sample_id]
            else:
                logging.getLogger('opal').warning(
                    "Skipping assessment of {} for sample {}. Make sure the SampleID of the gold standard and the profile are identical.\n".format(
                        label, sample_id))
                continue

            rank_to_taxid_to_percentage = load_data.get_rank_to_taxid_to_percentage(profile)

            unifrac, unifrac_cami, shannon, l1norm, binary_metrics, braycurtis, rank_to_sum, rank_to_ntaxa = \
                compute_metrics(sample_metadata, profile, gs_pf_profile,
                                profile, gs_pf_profile_cami,
                                gs_rank_to_taxid_to_percentage,
                                rank_to_taxid_to_percentage,
                                branch_length_fun, normalized_unifrac)
            rename_as_unfiltered = True if filter_tail_percentage else False
            pd_metrics = pd.concat([pd_metrics,
                                    reformat_pandas(sample_id, label, braycurtis, shannon, binary_metrics, l1norm,
                                                    unifrac, unifrac_cami, rank_to_sum, rank_to_ntaxa,
                                                    rename_as_unfiltered)], ignore_index=True)
            pd_sample_confusion = get_confusion_df(gs_rank_to_taxid_to_percentage, rank_to_taxid_to_percentage)
            pd_sample_confusion['sample'] = sample_id
            pd_sample_confusion['tool'] = label
            pd_confusion = pd.concat([pd_confusion, pd_sample_confusion], ignore_index=True)

            if filter_tail_percentage:
                rank_to_taxid_to_percentage_filtered = \
                    load_data.get_rank_to_taxid_to_percentage_filtered(rank_to_taxid_to_percentage,
                                                                       filter_tail_percentage)
                profile_filtered = [prediction for prediction in profile if
                                    prediction.taxid in rank_to_taxid_to_percentage_filtered[prediction.rank]]
                unifrac, unifrac_cami, shannon, l1norm, binary_metrics, braycurtis, rank_to_sum, rank_to_ntaxa = \
                    compute_metrics(sample_metadata, profile_filtered, gs_pf_profile,
                                    profile_filtered, gs_pf_profile_cami,
                                    gs_rank_to_taxid_to_percentage,
                                    rank_to_taxid_to_percentage_filtered,
                                    branch_length_fun, normalized_unifrac)
                pd_metrics = pd.concat([pd_metrics,
                                        reformat_pandas(sample_id, label, braycurtis, shannon, binary_metrics, l1norm,
                                                        unifrac, unifrac_cami, rank_to_sum, rank_to_ntaxa)],
                                       ignore_index=True)

            one_profile_assessed = True

    if not one_profile_assessed:
        logging.getLogger('opal').critical("No profile could be evaluated.")
        exit(1)

    return pd_metrics, pd_confusion
