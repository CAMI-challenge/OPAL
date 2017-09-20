#!/usr/bin/env python

import os
import errno
import argparse
import os.path
import l1norm as l1
import binary_metrics as bm
import unifrac_distance as uf
from utils import load_data
from utils import ProfilingTools as PF


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


def evaluate(gold_standard_file, profiles_files, labels, output_dir):
    # L1 Norm
    l1norm_list = []
    weighted_unifrac_list = []
    gs_rank_to_taxid_to_percentage = load_data.open_profile(gold_standard_file)
    gs_profile = PF.Profile(input_file_name=gold_standard_file)
    for profile_file, label in zip(profiles_files, labels):
        rank_to_taxid_to_percentage = load_data.open_profile(profile_file)
        pf_profile = PF.Profile(input_file_name=profile_file)
        l1norm_list.append(l1.compute_l1norm(gs_rank_to_taxid_to_percentage, rank_to_taxid_to_percentage))
        compute_binary_metrics(rank_to_taxid_to_percentage, gs_rank_to_taxid_to_percentage, os.path.join(output_dir, label))
        weighted_unifrac_list.append(uf.compute_unifrac(gs_profile, pf_profile))

    f = open(output_dir + "/l1_norm.tsv", 'w')
    l1.print_list_l1norm(l1norm_list, labels, f)
    f.close()

    f = open(output_dir + "/unifrac.tsv", 'w')
    uf.print_list_unifrac(weighted_unifrac_list, labels, f)
    f.close()


def main():
    parser = argparse.ArgumentParser(description="Compute all metrics for one or more profiles")
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
