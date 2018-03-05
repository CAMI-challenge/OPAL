#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
import itertools
from collections import defaultdict
from collections import OrderedDict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from mpl_toolkits.mplot3d import Axes3D
import braycurtis as bc
from utils import spider_plot_functions as spl
from utils import constants as c
from utils import load_data


def create_colors_list():
    colors_list = []
    for color in plt.cm.tab10(np.linspace(0, 1, 10))[:-1]:
        colors_list.append(tuple(color))
    colors_list.append("black")
    for color in plt.cm.Set2(np.linspace(0, 1, 8)):
        colors_list.append(tuple(color))
    for color in plt.cm.Set3(np.linspace(0, 1, 12)):
        colors_list.append(tuple(color))
    return colors_list


def do_scatter_plot(profile_values, gs_values, output_dir, label):
    fig, axs = plt.subplots(figsize=(6, 5))

    # force axes to be from 0 to 100%
    axs.set_xlim([0.0, 1.0])
    axs.set_ylim([0.0, 1.0])

    # turn on grid
    axs.grid(which='major', linestyle='-', linewidth='0.5', color='lightgrey')

    axs.scatter(profile_values, gs_values, alpha=0.5)
    axs.set_xlabel(label)
    axs.set_ylabel('Gold standard')
    axs.set_title(c.BRAY_CURTIS)

    fig.savefig(os.path.join(output_dir, 'beta_diversity_bc_' + label + '.pdf'), dpi=100, format='pdf', bbox_inches='tight')


def plot_beta_diversity(gs_samples_list, profiles_list_to_samples_list, sample_ids_list, labels, output_dir):
    if len(sample_ids_list) < 2:
        return

    gs_sampleid_to_rank_to_taxid_to_percentage = {}
    for sample in gs_samples_list:
        sample_id, sample_metadata, profile = sample
        gs_sampleid_to_rank_to_taxid_to_percentage[sample_id] = load_data.get_rank_to_taxid_to_percentage(profile, c.SPECIES)

    sample_ids_combinations = list(itertools.combinations(sample_ids_list, 2))
    gs_braycurtis = {}
    # Compute Bray-Curtis for all combinations of samples
    for sample1, sample2 in sample_ids_combinations:
        gs_braycurtis[(sample1, sample2)] = bc.braycurtis(gs_sampleid_to_rank_to_taxid_to_percentage[sample1], gs_sampleid_to_rank_to_taxid_to_percentage[sample2])[c.SPECIES]

    profile_sampleid_to_rank_to_taxid_to_percentage = {}
    for profile, label in zip(profiles_list_to_samples_list, labels):
        profile_sampleid_to_rank_to_taxid_to_percentage[label] = {}
        for sample in profile:
            sample_id, sample_metadata, profile = sample
            if not sample_id in gs_sampleid_to_rank_to_taxid_to_percentage:
                continue
            profile_sampleid_to_rank_to_taxid_to_percentage[label][sample_id] = load_data.get_rank_to_taxid_to_percentage(profile, c.SPECIES)

        profile_braycurtis = {}
        # Compute Bray-Curtis for all combinations of samples
        for sample1, sample2 in sample_ids_combinations:
            if sample1 in profile_sampleid_to_rank_to_taxid_to_percentage[label] and sample2 in profile_sampleid_to_rank_to_taxid_to_percentage[label]:
                profile_braycurtis[(sample1, sample2)] = bc.braycurtis(profile_sampleid_to_rank_to_taxid_to_percentage[label][sample1], profile_sampleid_to_rank_to_taxid_to_percentage[label][sample2])[c.SPECIES]

        gs_values = []
        profile_values = []
        for (sample1, sample2), value in profile_braycurtis.items():
            gs_values.append(gs_braycurtis[(sample1, sample2)])
            profile_values.append(profile_braycurtis[(sample1, sample2)])
        do_scatter_plot(profile_values, gs_values, output_dir, label)


def plot_shannon(rank_to_shannon_list, rank_to_shannon_gs, labels, output_dir, file_name):
    colors_list = create_colors_list()

    if rank_to_shannon_gs is None:
        # skip color of gold standard, if not present
        colors_list = colors_list[1:]
        ylabel = 'Shannon equitability (absolute diff. to gold standard)'
    else:
        labels = ['Gold standard'] + labels
        ylabel = 'Shannon equitability'

    fig, axs = plt.subplots(figsize=(6, 5))

    # force axis to be from 0 to 100%
    axs.set_ylim([0.0, 1.0])

    # turn on grid
    axs.grid(which='major', linestyle='-', linewidth='0.5', color='lightgrey')

    for i, rank_to_shannon in enumerate(rank_to_shannon_gs + rank_to_shannon_list if rank_to_shannon_gs is not None else rank_to_shannon_list):
        x = []
        y = []
        for j, rank in enumerate(c.ALL_RANKS, start=1):
            if rank in rank_to_shannon:
                x.append(j)
                y.append(rank_to_shannon[rank])
        axs.plot(x, y, color=colors_list[i], marker='o', markersize=10)

    axs.set_xticklabels([''] + c.ALL_RANKS)
    plt.setp(axs.get_xticklabels(), fontsize=9, rotation=25)

    axs.set_ylabel(ylabel)

    lgd = plt.legend(labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., handlelength=0, frameon=False)
    fig.savefig(os.path.join(output_dir, file_name + '.pdf'), dpi=100, format='pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, file_name + '.png'), dpi=100, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')


def spider_plot(metrics, labels, rank_to_metric_to_toolvalues, output_dir, file_name, colors, grid_points=None, fill=False):
    N = len(labels)
    if N < 3:
        return
    theta = spl.radar_factory(N, frame='polygon')
    fig, axes = plt.subplots(figsize=(9, 9), nrows=2, ncols=3, subplot_kw=dict(projection='radar'))
    fig.subplots_adjust(wspace=0.35, hspace=0.05, top=0.87, bottom=0.3)

    for ax, rank in zip(axes.flatten(), c.PHYLUM_SPECIES):
        if grid_points:
            ax.set_rgrids(grid_points, fontsize='xx-small')
        else:
            ax.set_rgrids([0.2, 0.4, 0.6, 0.8], ('', '', '', ''))  # get rid of the labels of the grid points
        ax.set_title(rank, weight='bold', size='medium', position=(0.5, 1.1),
                     horizontalalignment='center', verticalalignment='center')

        # select only metrics in metrics list
        metrics_subdict = OrderedDict((metric, rank_to_metric_to_toolvalues[rank][metric]) for metric in metrics)
        it = 1
        metric_to_toolindex = []
        for d, color in zip(metrics_subdict.values(), colors):
            # store index of tools without a value for the current metric
            metric_to_toolindex.append([i for i, x in enumerate(d) if x is None or np.isnan(x)])
            d = [0 if x is None or np.isnan(x) else x for x in d]

            ax.plot(theta, d, '--', color=color, linewidth=3, dashes=(it, 1))
            if fill:
                ax.fill(theta, d, facecolor=color, alpha=0.25)
            it += 1
        ax.set_varlabels(labels)

        ax.set_ylim([0.0, 1.0])
        ax.set_xlim([0.0, 1.0])

        # color red label of tools without a value for at least one metric
        xticklabels = ax.get_xticklabels()
        for metric in metric_to_toolindex:
            for toolindex in metric:
                xticklabels[toolindex].set_color([1, 0, 0])

    ax = axes[0, 0]
    ax.legend(metrics, loc=(1.68 - 0.353 * len(metrics), 1.3), labelspacing=0.1, fontsize='small', ncol=len(metrics))
    fig.savefig(os.path.join(output_dir, file_name + '.pdf'), dpi=100, format='pdf', bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, file_name + '.png'), dpi=100, format='png', bbox_inches='tight')


def plot_braycurtis_l1norm(braycurtis_list, l1norm_list, labels, output_dir):
    colors_list = create_colors_list()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    reversed_ranks = c.ALL_RANKS.copy()
    reversed_ranks.reverse()
    rank_index = []

    color_index = 0
    for braycurtis, l1norm in zip(braycurtis_list, l1norm_list):
        for i, rank in enumerate(reversed_ranks, start=1):
            if rank in braycurtis and rank in l1norm:
                ax.scatter(braycurtis[rank], l1norm[rank], i, color=colors_list[color_index])
                rank_index.append(i)
        color_index += 1

    ax.set_xlabel('Bray-Curtis')
    ax.set_ylabel('L1 norm')

    rank_index = list(sorted(set(rank_index)))
    zlabels = [reversed_ranks[i - 1] for i in rank_index]
    ax.set_zticklabels(zlabels)

    legend_handles = []
    for i, label in enumerate(labels):
        line = mlines.Line2D([], [], color='white', markerfacecolor=colors_list[i], marker='o', markersize=8, label=label)
        legend_handles.append(line)

    plt.legend(handles=legend_handles, loc=2)

    fig.savefig(os.path.join(output_dir, 'braycurtis_l1norm.pdf'), dpi=100, format='pdf', bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, 'braycurtis_l1norm.png'), dpi=100, format='png', bbox_inches='tight')


def plot_all(pd_metrics, labels, output_dir):
    metrics = [c.UNIFRAC, c.L1NORM, c.RECALL, c.PRECISION, c.FP]
    rank_to_metric_to_toolvalues = defaultdict(lambda : defaultdict(list))

    pd_copy = pd_metrics.copy()
    pd_copy.loc[pd_copy[pd.isnull(pd_copy['rank'])].index, 'rank'] = 'rank independent'
    # transform table
    pd_grouped = pd_copy.pivot_table(index=['rank', 'tool', 'sample'], columns='metric', values='value')

    # compute mean over samples and collapse rows
    cols = pd_grouped.columns
    pd_grouped = pd_grouped.groupby(['rank', 'tool'], sort=False)[cols].mean()

    # store maximum precision per rank
    rank_to_max_fp = defaultdict()
    for rank in c.PHYLUM_SPECIES:
        pd_rank = pd_grouped.loc[(pd_grouped.index.get_level_values('rank') == rank)]
        rank_to_max_fp[rank] = pd_rank[c.FP].max()

    tool_to_rank_to_metric_to_value = defaultdict(lambda: defaultdict(dict))
    for (rank, tool), g in pd_grouped.groupby(['rank', 'tool']):
        if tool not in labels:
            continue
        if rank == 'rank independent':
            tool_to_rank_to_metric_to_value[tool][rank][c.UNIFRAC] = g[c.UNIFRAC].values[0] /16
        elif rank in c.PHYLUM_SPECIES:
            tool_to_rank_to_metric_to_value[tool][rank][c.L1NORM] = g[c.L1NORM].values[0] / 2.0 if len(g[c.L1NORM].values) > 0 else None
            tool_to_rank_to_metric_to_value[tool][rank][c.RECALL] = g[c.RECALL].values[0] if len(g[c.RECALL].values) > 0 else None
            tool_to_rank_to_metric_to_value[tool][rank][c.PRECISION] = g[c.PRECISION].values[0] if len(g[c.PRECISION].values) > 0 else None
            if rank_to_max_fp[rank] > 0:
                tool_to_rank_to_metric_to_value[tool][rank][c.FP] = g[c.FP].values[0] / rank_to_max_fp[rank] if len(g[c.FP].values) else None
            else:
                tool_to_rank_to_metric_to_value[tool][rank][c.FP] = g[c.FP].values[0] if len(g[c.FP].values) else None

    present_labels = []
    for label in labels:
        if label not in tool_to_rank_to_metric_to_value:
            continue
        else:
            present_labels.append(label)
        for rank in c.PHYLUM_SPECIES:
            for metric in metrics:
                if metric in tool_to_rank_to_metric_to_value[label][rank]:
                    rank_to_metric_to_toolvalues[rank][metric].append(tool_to_rank_to_metric_to_value[label][rank][metric])
            rank_to_metric_to_toolvalues[rank][c.UNIFRAC].append(tool_to_rank_to_metric_to_value[label]['rank independent'][c.UNIFRAC])

    spider_plot(metrics,
                present_labels,
                rank_to_metric_to_toolvalues,
                output_dir,
                'spider_plot',
                ['b', 'g', 'r', 'k', 'm'])

    spider_plot([c.RECALL, c.PRECISION],
                present_labels,
                rank_to_metric_to_toolvalues,
                output_dir,
                'spider_plot_recall_precision',
                ['r', 'k'],
                grid_points=[0.2, 0.4, 0.6, 0.8, 1.0],
                fill=True)

    # compute average shannon for gold standard
    pd_shannon_equit = pd_metrics[pd_metrics['metric'] == c.SHANNON_EQUIT]
    table1 = pd_shannon_equit[pd_shannon_equit['tool'] == c.GS][['rank', 'value']]
    rank_to_shannon_gs = table1.groupby('rank').mean().T.to_dict('records')

    tool_to_rank_to_shannon = defaultdict(dict)
    for (rank, tool), g in pd_grouped.groupby(['rank', 'tool']):
        if tool == c.GS or rank == 'rank independent':
            continue
        tool_to_rank_to_shannon[tool][rank] = g[c.SHANNON_EQUIT].values[0] if len(g[c.SHANNON_EQUIT].values) > 0 else None

    # compute differences between each tool and gold standard
    tool_to_rank_to_shannon_difference = defaultdict(dict)

    for tool in tool_to_rank_to_shannon.keys():
        for rank in tool_to_rank_to_shannon[tool].keys():
            if rank in rank_to_shannon_gs[0]:
                tool_to_rank_to_shannon_difference[tool][rank] = abs(rank_to_shannon_gs[0][rank] - tool_to_rank_to_shannon[tool][rank])

    # convert to list of dictionaries
    tool_to_rank_to_shannon = [tool_to_rank_to_shannon[label] for label in labels]
    tool_to_rank_to_shannon_difference = [tool_to_rank_to_shannon_difference[label] for label in labels]

    plot_shannon(tool_to_rank_to_shannon, rank_to_shannon_gs, labels, output_dir, 'plot_shannon')
    plot_shannon(tool_to_rank_to_shannon_difference, None, labels, output_dir, 'plot_shannon_diff')

    # pl.plot_braycurtis_l1norm(braycurtis_list, l1norm_list, labels, output_dir)
