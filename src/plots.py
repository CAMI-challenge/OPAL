#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
import itertools
import math
import copy
from collections import defaultdict
from collections import OrderedDict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.lines import Line2D
from src import braycurtis as bc
from src.utils import spider_plot_functions as spl
from src.utils import constants as c
from src.utils import load_data
import scipy.special


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


def plot_time_memory(time, memory, labels, output_dir):
    time_label = 'Time (hours)'
    memory_label = 'Memory (GB)'

    df = pd.DataFrame(index=labels)
    if time:
        df[time_label] = time
    if memory:
        df[memory_label] = memory

    if time:
        df.sort_values(by=[time_label], inplace=True)
    else:
        df.sort_values(by=[memory_label], inplace=True)

    label2 = ''
    if time:
        label1 = time_label
        if memory:
            label2 = memory_label
    else:
        label1 = memory_label

    if label2:
        ax = df.plot(kind='bar', secondary_y=label2, mark_right=False, zorder=20)
    else:
        ax = df.plot(kind='bar', mark_right=False, zorder=20)
    ax.grid(which='major', linestyle='-', linewidth='0.5', color='lightgrey', zorder=0)
    ax.set_ylabel(label1)
    plt.setp(ax.get_xticklabels(), fontsize=9, rotation=15)

    if label2:
        ax2 = ax.twinx()
        ax2.set_ylabel(memory_label, labelpad=23)
        ax2.yaxis.set_ticks([])
    fig_name = 'time_memory'
    plt.savefig(os.path.join(output_dir, fig_name + '.pdf'), dpi=100, format='pdf', bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, fig_name + '.png'), dpi=100, format='png', bbox_inches='tight')
    return [fig_name]


def create_legend_rarefaction(output_dir):
    colors_iter = iter(create_colors_list())
    circles = []
    for x in c.ALL_RANKS:
        nextcolor = next(colors_iter)
        circles.append(Line2D([], [], markeredgewidth=2.0, linestyle="None", marker=0, markersize=15, markeredgecolor=nextcolor, markerfacecolor=nextcolor))

    fig = plt.figure(figsize=(0.5, 0.5))
    fig.legend(circles, c.ALL_RANKS, loc='center', frameon=False, ncol=8, handletextpad=-0.5, columnspacing=2.0)
    fig.savefig(os.path.join(output_dir, 'gold_standard', 'legend.eps'), dpi=100, format='eps', bbox_inches='tight')
    plt.close(fig)


def plot_rarefaction_curves(gs_samples_list, output_dir, log_scale=False):
    colors_list = create_colors_list()
    fig, axs = plt.subplots(figsize=(6, 5))

    # turn on grid
    axs.grid(which='major', linestyle='-', linewidth='0.5', color='lightgrey')

    axs.xaxis.set_major_locator(MaxNLocator(integer=True))

    axs.set_xlabel('Number of samples')
    if log_scale:
        axs.set_ylabel('log$_{10}$(number of taxa)')
    else:
        axs.set_ylabel('Number of taxa')

    x = list(range(1, len(gs_samples_list) + 1))
    y_rank_to_rar = OrderedDict((rank, []) for rank in c.ALL_RANKS)
    y_rank_to_acc = OrderedDict((rank, []) for rank in c.ALL_RANKS)

    for i, rank in enumerate(c.ALL_RANKS):
        tax_id_to_num_occurrences = {}

        # compute accumulation curves
        num_ids = 0
        for sample in gs_samples_list:
            sample_id, sample_metadata, profile = sample
            for prediction in profile:
                if prediction.rank == rank and prediction.percentage > .0:
                    if prediction.taxid in tax_id_to_num_occurrences:
                        tax_id_to_num_occurrences[prediction.taxid] += 1
                    else:
                        tax_id_to_num_occurrences[prediction.taxid] = 1
            num_ids = len(tax_id_to_num_occurrences)
            if log_scale:
                y_rank_to_rar[rank].append(math.log10(num_ids) if num_ids > 0 else 0)
            else:
                y_rank_to_rar[rank].append(num_ids)

        # compute rarefaction curves
        num_samples = len(gs_samples_list)
        for j, sample in enumerate(gs_samples_list, 1):
            sum_over_taxa = 0
            for tax_id in tax_id_to_num_occurrences:
                sum_over_taxa += scipy.special.binom(num_samples - tax_id_to_num_occurrences[tax_id], j)
            if log_scale:
                value = num_ids - (scipy.special.binom(num_samples, j) ** (-1) * sum_over_taxa)
                y_rank_to_acc[rank].append(math.log10(value) if value > 0 else 0)
            else:
                y_rank_to_acc[rank].append(num_ids - (scipy.special.binom(num_samples, j) ** (-1) * sum_over_taxa))

    for i, rank in enumerate(c.ALL_RANKS):
        axs.plot(x, y_rank_to_acc[rank], color=colors_list[i])
    for i, rank in enumerate(c.ALL_RANKS):
        axs.plot(x, y_rank_to_rar[rank], linestyle=':', color=colors_list[i])

    if log_scale:
        file_name = 'rarefaction_curves_log_scale'
    else:
        file_name = 'rarefaction_curves'

    # fig.savefig(os.path.join(output_dir, 'gold_standard', file_name + '_wolegend.pdf'), dpi=100, format='pdf', bbox_inches='tight')

    lgd = plt.legend(c.ALL_RANKS, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., handlelength=2, frameon=False)
    fig.savefig(os.path.join(output_dir, 'gold_standard', file_name + '.pdf'), dpi=100, format='pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, 'gold_standard', file_name + '.png'), dpi=100, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)
    return [file_name]


def plot_samples_hist(gs_samples_list, sample_ids_list, output_dir):
    plots_list = []
    for rank in c.ALL_RANKS:
        df_gs = pd.DataFrame()
        for sample in gs_samples_list:
            sample_id, sample_metadata, profile = sample
            gs_sampleid_to_rank_to_taxid_to_percentage = {sample_id: load_data.get_rank_to_taxid_to_percentage(profile, rank)}
            if len(gs_sampleid_to_rank_to_taxid_to_percentage[sample_id]) == 0:
                continue

            # add taxa names
            tax_id_to_name = load_data.get_taxa_names(profile)
            gs_sampleid_to_rank_to_taxidandname_to_percentage = copy.deepcopy(gs_sampleid_to_rank_to_taxid_to_percentage)
            for tax_id, percentage in gs_sampleid_to_rank_to_taxid_to_percentage[sample_id][rank].items():
                gs_sampleid_to_rank_to_taxidandname_to_percentage[sample_id][rank][tax_id_to_name[tax_id] + ' ' + tax_id] = percentage
                del gs_sampleid_to_rank_to_taxidandname_to_percentage[sample_id][rank][tax_id]

            df_gs_sample = pd.DataFrame.from_dict(gs_sampleid_to_rank_to_taxidandname_to_percentage[sample_id][rank], orient='index')
            df_gs_sample.index.name = 'taxid'
            df_gs_sample['Sample'] = sample_id
            df_gs_sample.rename(columns={0: 'value'}, inplace=True)
            df_gs_sample.set_index('Sample', append=True, inplace=True)
            df_gs = pd.concat([df_gs, df_gs_sample])
        if df_gs.empty:
            continue

        df_gs = df_gs.pivot_table(index=['taxid'], columns='Sample', values='value').T

        df_gs = df_gs.loc[sample_ids_list]

        fig, axs = plt.subplots(figsize=(6, 5))
        axs.tick_params(axis='x',
                        which='both',
                        bottom='off',
                        top='off',
                        labelbottom='off')
        axs.set_ylabel('Proportions [%]')
        axs.set_title('Rank: ' + rank)
        fig_name = os.path.join('gold_standard', rank)
        fig = df_gs.plot(kind='bar', ax=axs, stacked=True, legend=False, colormap=plt.get_cmap('gist_rainbow')).get_figure()
        fig.savefig(os.path.join(output_dir, fig_name + '.pdf'), dpi=100, format='pdf', bbox_inches='tight')
        fig.savefig(os.path.join(output_dir, fig_name + '.png'), dpi=100, format='png', bbox_inches='tight')

        legend, labels = axs.get_legend_handles_labels()
        fig2, axs2 = plt.subplots(figsize=(1, 1))
        axs2.xaxis.set_visible(False)
        axs2.yaxis.set_visible(False)
        for v in axs2.spines.values():
            v.set_visible(False)

        ncol = int(math.ceil(df_gs.shape[1]/35))
        axs2.legend(legend, labels, loc='center', frameon=True, ncol=ncol, handletextpad=0, columnspacing=1.0, borderaxespad=0., fontsize=9)
        plt.tight_layout()
        fig2.savefig(os.path.join(output_dir, fig_name + '_legend.png'), dpi=100, format='png', bbox_inches='tight')

        plt.close(fig)
        plt.close(fig2)
        plots_list.append(fig_name)
    return plots_list


def do_scatter_plot(profile_values, gs_values, output_dir, rank, label):
    fig, axs = plt.subplots(figsize=(6, 5))

    # force axes to be from 0 to 100%
    axs.set_xlim([0.0, 1.0])
    axs.set_ylim([0.0, 1.0])

    # turn on grid
    axs.grid(which='major', linestyle='-', linewidth='0.5', color='lightgrey')

    axs.scatter(profile_values, gs_values, alpha=0.5)
    axs.set_xlabel(label)
    axs.set_ylabel('Gold standard')
    axs.set_title(c.BRAY_CURTIS + ' - Rank: ' + rank)

    fig_name = os.path.join("by_tool", label, 'beta_diversity_bc_' + rank)
    fig.savefig(os.path.join(output_dir, fig_name + '.pdf'), dpi=100, format='pdf', bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, fig_name + '.png'), dpi=100, format='png', bbox_inches='tight')
    plt.close(fig)
    return [fig_name]


def plot_beta_diversity(gs_samples_list, profiles_list_to_samples_list, sample_ids_list, labels, output_dir):
    if len(sample_ids_list) < 2:
        return []

    gs_sampleid_to_rank_to_taxid_to_percentage = {}
    for sample in gs_samples_list:
        sample_id, sample_metadata, profile = sample
        gs_sampleid_to_rank_to_taxid_to_percentage[sample_id] = load_data.get_rank_to_taxid_to_percentage(profile)

    # Gold standard
    sample_ids_combinations = list(itertools.combinations(sample_ids_list, 2))
    gs_sample12_to_rank_to_braycurtis = {}
    # Compute Bray-Curtis for all combinations of samples
    for sample1, sample2 in sample_ids_combinations:
        gs_sample12_to_rank_to_braycurtis[(sample1, sample2)] = bc.braycurtis(gs_sampleid_to_rank_to_taxid_to_percentage[sample1], gs_sampleid_to_rank_to_taxid_to_percentage[sample2])
    # Sort by rank
    gs_rank_to_sample12_to_braycurtis = defaultdict(dict)
    for sample12, rank_to_braycurtis in gs_sample12_to_rank_to_braycurtis.items():
        for rank, value in rank_to_braycurtis.items():
            gs_rank_to_sample12_to_braycurtis[rank][sample12] = value

    plots_list = []

    # Assessed profiles
    profile_sampleid_to_rank_to_taxid_to_percentage = {}
    for profile, label in zip(profiles_list_to_samples_list, labels):
        profile_sampleid_to_rank_to_taxid_to_percentage[label] = {}
        for sample in profile:
            sample_id, sample_metadata, profile = sample
            if not sample_id in gs_sampleid_to_rank_to_taxid_to_percentage:
                continue
            profile_sampleid_to_rank_to_taxid_to_percentage[label][sample_id] = load_data.get_rank_to_taxid_to_percentage(profile)

        profile_sample12_to_rank_to_braycurtis = {}
        # Compute Bray-Curtis for all combinations of samples
        for sample1, sample2 in sample_ids_combinations:
            if sample1 in profile_sampleid_to_rank_to_taxid_to_percentage[label] and sample2 in profile_sampleid_to_rank_to_taxid_to_percentage[label]:
                profile_sample12_to_rank_to_braycurtis[(sample1, sample2)] = bc.braycurtis(profile_sampleid_to_rank_to_taxid_to_percentage[label][sample1], profile_sampleid_to_rank_to_taxid_to_percentage[label][sample2])
        # Sort by rank
        profile_rank_to_sample12_to_braycurtis = defaultdict(dict)
        for sample12, rank_to_braycurtis in profile_sample12_to_rank_to_braycurtis.items():
            for rank, value in rank_to_braycurtis.items():
                profile_rank_to_sample12_to_braycurtis[rank][sample12] = value

        for rank in c.ALL_RANKS:
            gs_values = []
            profile_values = []
            if rank in gs_rank_to_sample12_to_braycurtis and rank in profile_rank_to_sample12_to_braycurtis:
                for (sample1, sample2), value in profile_rank_to_sample12_to_braycurtis[rank].items():
                    gs_values.append(gs_rank_to_sample12_to_braycurtis[rank][(sample1, sample2)])
                    profile_values.append(profile_rank_to_sample12_to_braycurtis[rank][(sample1, sample2)])
                plots_list += do_scatter_plot(profile_values, gs_values, output_dir, rank, label)
    return plots_list


def create_legend_shannon(labels, output_dir):
    colors_iter = iter(create_colors_list())
    labels = ['Gold standard'] + labels
    circles = []
    for x in c.ALL_RANKS:
        nextcolor = next(colors_iter)
        circles.append(Line2D([], [], markeredgewidth=2.0, linestyle="None", marker="o", markersize=6, markeredgecolor=nextcolor, markerfacecolor=nextcolor))

    fig = plt.figure(figsize=(0.5, 0.5))
    fig.legend(circles, labels, loc='center', frameon=False, ncol=8, handletextpad=0, columnspacing=1.0)
    fig.savefig(os.path.join(output_dir, 'legend_shannon.eps'), dpi=100, format='eps', bbox_inches='tight')
    plt.close(fig)


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

    # fig.savefig(os.path.join(output_dir, file_name + 'wolegend.pdf'), dpi=100, format='pdf', bbox_inches='tight')

    lgd = plt.legend(labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., handlelength=0, frameon=False)
    fig.savefig(os.path.join(output_dir, file_name + '.pdf'), dpi=100, format='pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, file_name + '.png'), dpi=100, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close(fig)


def spider_plot(metrics, labels, rank_to_metric_to_toolvalues, output_dir, file_name, colors, grid_points=None, fill=False, absolute=False):
    N = len(labels)
    if N < 3:
        return []
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

        if absolute:
            metric_suffix = 'absolute'
        else:
            metric_suffix = ''
        # select only metrics in metrics list
        metrics_subdict = OrderedDict((metric, rank_to_metric_to_toolvalues[rank][metric + metric_suffix]) for metric in metrics)
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
    ax.legend(metrics, loc=(1.6 - 0.353 * len(metrics), 1.3), labelspacing=0.1, fontsize='small', ncol=len(metrics))
    fig.savefig(os.path.join(output_dir, file_name + '.pdf'), dpi=100, format='pdf', bbox_inches='tight')
    fig.savefig(os.path.join(output_dir, file_name + '.png'), dpi=100, format='png', bbox_inches='tight')
    plt.close(fig)
    return [file_name]


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
    plt.close(fig)


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

    # store maximum values per rank
    rank_to_max_fp = defaultdict()
    rank_to_max_unifrac = defaultdict()
    rank_to_max_l1norm = defaultdict()
    rank_to_max_recall = defaultdict()
    rank_to_max_precision = defaultdict()
    for rank in c.PHYLUM_SPECIES:
        pd_rank = pd_grouped.loc[(pd_grouped.index.get_level_values('rank') == rank) & (pd_grouped.index.get_level_values('tool') != c.GS)]
        rank_to_max_fp[rank] = pd_rank[c.FP].max()
        rank_to_max_l1norm[rank] = pd_rank[c.L1NORM].max()
        rank_to_max_recall[rank] = pd_rank[c.RECALL].max()
        rank_to_max_precision[rank] = pd_rank[c.PRECISION].max()
    pd_rank = pd_grouped.loc[(pd_grouped.index.get_level_values('rank') == 'rank independent') & (pd_grouped.index.get_level_values('tool') != c.GS)]
    rank_to_max_unifrac['rank independent'] = pd_rank[c.UNIFRAC].max()

    tool_to_rank_to_metric_to_value = defaultdict(lambda: defaultdict(dict))
    for (rank, tool), g in pd_grouped.groupby(['rank', 'tool']):
        if tool not in labels:
            continue
        if rank == 'rank independent':
            # relative values
            tool_to_rank_to_metric_to_value[tool][rank][c.UNIFRAC] = (g[c.UNIFRAC].values[0] / (rank_to_max_unifrac[rank]) if rank_to_max_unifrac[rank] > 0 else g[c.UNIFRAC].values[0])
        elif rank in c.PHYLUM_SPECIES:
            # absolute values
            tool_to_rank_to_metric_to_value[tool][rank][c.RECALL+'absolute'] = g[c.RECALL].values[0] if len(g[c.RECALL].values) > 0 else None
            tool_to_rank_to_metric_to_value[tool][rank][c.PRECISION+'absolute'] = g[c.PRECISION].values[0] if len(g[c.PRECISION].values) > 0 else None

            # relative values
            if rank_to_max_recall[rank] > 0:
                tool_to_rank_to_metric_to_value[tool][rank][c.RECALL] = (g[c.RECALL].values[0] / rank_to_max_recall[rank]) if len(g[c.RECALL].values) > 0 else None
            else:
                tool_to_rank_to_metric_to_value[tool][rank][c.RECALL] = g[c.RECALL].values[0] if len(g[c.RECALL].values) > 0 else None

            if rank_to_max_precision[rank] > 0:
                tool_to_rank_to_metric_to_value[tool][rank][c.PRECISION] = (g[c.PRECISION].values[0] / rank_to_max_precision[rank]) if len(g[c.PRECISION].values) > 0 else None
            else:
                tool_to_rank_to_metric_to_value[tool][rank][c.PRECISION] = g[c.PRECISION].values[0] if len(g[c.PRECISION].values) > 0 else None

            if rank_to_max_l1norm[rank] > 0:
                tool_to_rank_to_metric_to_value[tool][rank][c.L1NORM] = (g[c.L1NORM].values[0] / rank_to_max_l1norm[rank]) if len(g[c.L1NORM].values) > 0 else None
            else:
                tool_to_rank_to_metric_to_value[tool][rank][c.L1NORM] = g[c.L1NORM].values[0] if len(g[c.L1NORM].values) > 0 else None

            if rank_to_max_fp[rank] > 0:
                tool_to_rank_to_metric_to_value[tool][rank][c.FP] = (g[c.FP].values[0] / rank_to_max_fp[rank]) if len(g[c.FP].values) else None
            else:
                tool_to_rank_to_metric_to_value[tool][rank][c.FP] = g[c.FP].values[0] if len(g[c.FP].values) else None

    present_labels = []
    for label in labels:
        if label not in tool_to_rank_to_metric_to_value:
            continue
        else:
            present_labels.append(label)
        for rank in c.PHYLUM_SPECIES:
            for metric in metrics + [c.RECALL+'absolute', c.PRECISION+'absolute']:
                if metric in tool_to_rank_to_metric_to_value[label][rank]:
                    rank_to_metric_to_toolvalues[rank][metric].append(tool_to_rank_to_metric_to_value[label][rank][metric])
            rank_to_metric_to_toolvalues[rank][c.UNIFRAC].append(tool_to_rank_to_metric_to_value[label]['rank independent'][c.UNIFRAC])

    plots_list = spider_plot(metrics,
                             present_labels,
                             rank_to_metric_to_toolvalues,
                             output_dir,
                             'spider_plot',
                             ['b', 'g', 'r', 'k', 'm'])

    plots_list += spider_plot([c.RECALL, c.PRECISION],
                              present_labels,
                              rank_to_metric_to_toolvalues,
                              output_dir,
                              'spider_plot_recall_precision',
                              ['r', 'k'],
                              grid_points=[0.2, 0.4, 0.6, 0.8, 1.0],
                              fill=True,
                              absolute=True)

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
    return plots_list
