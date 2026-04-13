#!/usr/bin/env python

from cami_opal.utils import constants as c
import pandas as pd
import logging


def get_user_ranks_list(ranks):
    rank_high_low = [x.strip() for x in ranks.split(',')]
    if len(rank_high_low) != 2 or rank_high_low[0] not in c.RANK_GROUPS or rank_high_low[1] not in c.RANK_GROUPS:
        logging.getLogger('opal').warning('Invalid ranks provided with option --ranks. Default will be used.')
        return list(c.RANK_GROUPS.keys())[:-1]
    all_ranks = list(c.RANK_GROUPS.keys())
    index1 = all_ranks.index(rank_high_low[0])
    index2 = all_ranks.index(rank_high_low[1])
    if index1 < index2:
        return all_ranks[index1:index2 + 1]
    else:
        return all_ranks[index2:index1 + 1]


def highscore_table(metrics, ranks, collapse=False):
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
        Old default (CAMI 1): 'phylum', 'class', 'order', 'family', 'genus'
        Which ranks should be considered for rank dependent metrics.
        Here we decided to exclude e.g. species, because most profilers
        fail at that rank and we don't want to emphasize on this rank.
    Returns
    -------
    Pandas.DataFrame holding a high scoring table as in Figure 3c.
    """
    if ranks:
        useranks = get_user_ranks_list(ranks)
    else:
        useranks = list(dict.fromkeys(c.RANK_GROUPS.values()))[:-1] if collapse else list(c.RANK_GROUPS.keys())[:-1]

    pd_metrics = metrics.copy()
    pd_metrics.loc[pd_metrics[pd.isnull(pd_metrics['rank'])].index, 'rank'] = 'rank independent'

    sort_ascendingly = {c.L1NORM: True, c.UNIFRAC: True,
                        c.RECALL: False, c.PRECISION: False}

    # get rid of gold standard
    pd_metrics = pd_metrics[pd_metrics['tool'] != c.GS]

    # collecting rank scores
    posresults = []
    for (metric, sample, rank), g in pd_metrics.groupby(['metric', 'sample', 'rank']):
        if metric in sort_ascendingly:
            if ((rank in useranks) and (metric != c.UNIFRAC)) or ((rank == 'rank independent') and (metric == c.UNIFRAC)):
                res = g.groupby('tool').sum(numeric_only=True)
                res['position'] = res['value'].rank(method='min', ascending=sort_ascendingly[metric]) - 1
                res['metric'] = metric
                res['sample'] = sample
                res['rank'] = rank
                posresults.append(res)
    posresults = pd.concat(posresults)

    return posresults.groupby(['metric', 'tool'])['position'].sum().to_frame(), useranks

    # reformat like Figure 3c
    os = []
    for metric, g in posresults.groupby('metric'):
        highscore = g.groupby('tool')['position'].sum().sort_values()
        os.append(pd.DataFrame(["%s (%i)" % (idx, pos) for idx, pos in highscore.iteritems()], columns=[metric]))

    # return reformatted table
    return pd.concat(os, axis=1)  #.T.loc[['Recall', 'Precision', 'L1 norm error', 'Unweighted Unifrac error'],:]