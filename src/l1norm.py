#!/usr/bin/env python

import sys
import argparse
from src.utils import load_data
from src.utils import constants as c


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
        stream.write("\t".join([format(l1norm[c.SUPERKINGDOM], '.3f') if c.SUPERKINGDOM in l1norm else "na",
                      format(l1norm[c.PHYLUM], '.3f') if c.PHYLUM in l1norm else "na",
                      format(l1norm[c.CLASS], '.3f') if c.CLASS in l1norm else "na",
                      format(l1norm[c.ORDER], '.3f') if c.ORDER in l1norm else "na",
                      format(l1norm[c.FAMILY], '.3f') if c.FAMILY in l1norm else "na",
                      format(l1norm[c.GENUS], '.3f') if c.GENUS in l1norm else "na",
                      format(l1norm[c.SPECIES], '.3f') if c.SPECIES in l1norm else "na",
                      format(l1norm[c.STRAIN], '.3f') if c.STRAIN in l1norm else "na"]))
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