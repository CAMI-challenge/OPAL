#!/usr/bin/env python

import sys
import argparse
from utils import load_data


def compute_l1norm(f1_rank_to_taxid_to_percentage, f2_rank_to_taxid_to_percentage):
    l1norm = {}
    for rank in f1_rank_to_taxid_to_percentage:
        if rank not in f2_rank_to_taxid_to_percentage:
            continue
        l1norm[rank] = .0
        taxa_union = set(f1_rank_to_taxid_to_percentage[rank].keys()).union(f2_rank_to_taxid_to_percentage[rank].keys())
        for taxid in taxa_union:
            percentage1 = f1_rank_to_taxid_to_percentage[rank][taxid] if taxid in f1_rank_to_taxid_to_percentage[rank] else .0
            percentage2 = f2_rank_to_taxid_to_percentage[rank][taxid] if taxid in f2_rank_to_taxid_to_percentage[rank] else .0
            l1norm[rank] += abs(percentage1 - percentage2)
        l1norm[rank] = l1norm[rank] / 100.0
    return l1norm


def print_l1norm(l1norm, stream=sys.stdout):
    for rank in l1norm:
        stream.write("{}\t{}\n".format(rank, l1norm[rank]))


def print_list_l1norm(l1norm_list, labels, stream=sys.stdout):
    labels_iterator = iter(labels)
    stream.write("tool\tsuperkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tstrain\n")
    for l1norm in l1norm_list:
        stream.write("{}\t".format(next(labels_iterator)))
        stream.write("\t".join([format(l1norm["superkingdom"], '.3f') if "superkingdom" in l1norm else "na",
                      format(l1norm["phylum"], '.3f') if "phylum" in l1norm else "na",
                      format(l1norm["class"], '.3f') if "class" in l1norm else "na",
                      format(l1norm["order"], '.3f') if "order" in l1norm else "na",
                      format(l1norm["family"], '.3f') if "family" in l1norm else "na",
                      format(l1norm["genus"], '.3f') if "genus" in l1norm else "na",
                      format(l1norm["species"], '.3f') if "species" in l1norm else "na",
                      format(l1norm["strain"], '.3f') if "strain" in l1norm else "na"]))
        stream.write("\n")


def main():
    parser = argparse.ArgumentParser(description="Compute L1 Norm")
    parser.add_argument('-1', '--file1', help="Profile 1", required=True)
    parser.add_argument('-2', '--file2', help="Profile 2", required=True)
    args = parser.parse_args()
    f1_rank_to_taxid_to_percentage = load_data.open_profile(args.file1)
    f2_rank_to_taxid_to_percentage = load_data.open_profile(args.file2)
    l1norm = compute_l1norm(f1_rank_to_taxid_to_percentage, f2_rank_to_taxid_to_percentage)
    print_l1norm(l1norm)


if __name__ == "__main__":
    main()