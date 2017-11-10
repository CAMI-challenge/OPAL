#!/usr/bin/env python

import numpy as np
from collections import OrderedDict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from mpl_toolkits.mplot3d import Axes3D
from utils import spider_plot_functions as spl
from utils import constants as c


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


def plot_shannon(rank_to_shannon_list, rank_to_shannon_gs, labels, output_dir):
    colors_list = create_colors_list()

    fig, axs = plt.subplots(figsize=(6, 5))

    # force axis to be from 0 to 100%
    axs.set_ylim([0.0, 1.0])

    # turn on grid
    axs.grid(which='major', linestyle='-', linewidth='0.5', color='lightgrey')

    for i, rank_to_shannon in enumerate([rank_to_shannon_gs] + rank_to_shannon_list):
        x = []
        y = []
        for j, rank in enumerate(c.ALL_RANKS, start=1):
            if rank in rank_to_shannon:
                x.append(j)
                y.append(rank_to_shannon[rank].equitability)
        axs.plot(x, y, color=colors_list[i], marker='o', markersize=10)

    axs.set_xticklabels([''] + c.ALL_RANKS)
    plt.setp(axs.get_xticklabels(), fontsize=9, rotation=25)

    axs.set_ylabel('Shannon equitability')

    lgd = plt.legend(['Gold standard'] + labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., handlelength=0, frameon=False)
    fig.savefig(output_dir + '/plot_shannon.pdf', dpi=100, format='pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    fig.savefig(output_dir + '/plot_shannon.png', dpi=100, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')


def plot(metrics, labels, rank_to_metric_to_toolvalues, output_dir, file_name, colors, grid_points=None, fill=False):
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
            metric_to_toolindex.append([i for i, x in enumerate(d) if x is None])
            d = [0 if x is None else x for x in d]

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
    ax.legend(metrics, loc=(2.067 - 0.353 * len(metrics), 1.3), labelspacing=0.1, fontsize='small', ncol=len(metrics))
    fig.savefig(output_dir + '/' + file_name + '.pdf', dpi=100, format='pdf', bbox_inches='tight')
    fig.savefig(output_dir + '/' + file_name + '.png', dpi=100, format='png', bbox_inches='tight')


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

    fig.savefig(output_dir + '/braycurtis_l1norm.pdf', dpi=100, format='pdf', bbox_inches='tight')
    fig.savefig(output_dir + '/braycurtis_l1norm.png', dpi=100, format='png', bbox_inches='tight')
