#!/usr/bin/env python

import os
import errno
import argparse
import os.path
from collections import defaultdict
from collections import OrderedDict
import l1norm as l1
import binary_metrics as bm
import unifrac_distance as uf
from utils import plot_functions as pl
from utils import load_data
from utils import ProfilingTools as PF
from utils import constants as c
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


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


def compute_binary_metrics(query_profile, query_truth, path):
    all_metrics = bm.compute_tree_metrics(query_profile, query_truth)
    bm.print_all_metrics(all_metrics, path)
    return all_metrics


def plot(metrics, labels, rank_to_metric_to_toolvalues, output_dir, file_name, colors, grid_points=None, fill=False):
    N = len(labels)
    if N < 3:
        return
    theta = pl.radar_factory(N, frame='polygon')
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
    fig.savefig(output_dir + '/' + file_name, dpi=100, format='pdf', bbox_inches='tight')


def evaluate(gold_standard_file, profiles_files, labels, output_dir):
    l1norm_list = []
    binary_metrics_list = []
    weighted_unifrac_list = []
    gs_rank_to_taxid_to_percentage = load_data.open_profile(gold_standard_file)
    gs_profile = PF.Profile(input_file_name=gold_standard_file)
    metrics = [c.UNIFRAC, c.L1NORM, c.RECALL, c.PRECISION, c.FP]

    for profile_file, label in zip(profiles_files, labels):
        rank_to_taxid_to_percentage = load_data.open_profile(profile_file)
        pf_profile = PF.Profile(input_file_name=profile_file)

        # L1 Norm
        l1norm = l1.compute_l1norm(gs_rank_to_taxid_to_percentage, rank_to_taxid_to_percentage)
        l1norm_list.append(l1norm)

        # Binary metrics
        binary_metrics = compute_binary_metrics(rank_to_taxid_to_percentage, gs_rank_to_taxid_to_percentage, os.path.join(output_dir, label))
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
            rank_to_metric_to_toolvalues[rank][c.FP].append(binary_metrics[rank].fp / max_fp if rank in binary_metrics else None)

    plot(metrics,
         labels,
         rank_to_metric_to_toolvalues,
         output_dir,
         'plot.pdf',
         ['b', 'g', 'r', 'k', 'm'])

    plot([c.RECALL, c.PRECISION],
         labels,
         rank_to_metric_to_toolvalues,
         output_dir,
         'plot_recall_precision.pdf',
         ['r', 'k'],
         grid_points=[0.2, 0.4, 0.6, 0.8, 1.0],
         fill=True)

    f = open(output_dir + "/l1_norm.tsv", 'w')
    l1.print_list_l1norm(l1norm_list, labels, f)
    f.close()

    f = open(output_dir + "/unifrac.tsv", 'w')
    uf.print_list_unifrac(weighted_unifrac_list, labels, f)
    f.close()


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
    evaluate(args.gold_standard_file, args.profiles_files, labels, output_dir)


if __name__ == "__main__":
    main()
