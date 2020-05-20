# This is an example of how to use the script
from src.utils import ProfilingTools as PF
from src.utils import EMDUnifrac as EMDU
import numpy as np
import argparse
import copy
import sys
import os


def compute_unifrac(pf1, pf2):
    P1 = copy.deepcopy(pf1)
    P2 = copy.deepcopy(pf2)

    # weighted:
    (Tint, lint, nodes_in_order, nodes_to_index, P, Q) = P1.make_unifrac_input_and_normalize(P2)
    (weighted, _) = EMDU.EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
    ##################################################################
    # attempt 1:
    # compute normalizing factor: divide by worst possible
    # this is the theoretical max weighted unifrac. If the gold standard doesn't have abundance sum to 1 at the lowest
    # rank, the maximum weighted unifrac may be smaller # TODO
    # this also assumes that pf1 is the gold standard  # FIXME: but this is how Fernando is calling it, so should be ok
    gs_tax_path = P1.sample_metadata['RANKS'].split('|')  # effectively gets the max depth of the gs tree
    # max weighted unifrac:
    # 2->(all weight at tips whose LCA is the root of the taxonomic tree (so mass must move all the way up, then down)
    # 100->(convert to percentage)
    # tax_path_to_branch_len->(use the branch lengths defined in EMDUnifrac)
    max_weighted_unifrac = 2 * 100 * P1.tax_path_to_branch_len(gs_tax_path, P1.branch_len_func, P1.root_len)
    # TODO: Note: this makes the weighted unifrac a bit uninformative, since no reasonable profilier will get anywhere near the absolute worst weighted unifrac
    weighted = weighted / float(max_weighted_unifrac)
    ##################################################################
    # attempt 2:
    # ||WP - WQ|| / ||WP||
    # Basically: relative weighted unifrac error.
    #(gs_weighted_unifrac, _) = EMDU.EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, np.zeros(P.shape))
    #(raw_weighted_unifrac, _) = EMDU.EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
    #weighted = raw_weighted_unifrac / float(gs_weighted_unifrac)
    # TODO: Note: This really doesn't change the results at all, the same basic pattern is still there. revert to attempt 1

    # Unweighted
    #P1 = copy.deepcopy(pf1)
    #P2 = copy.deepcopy(pf2)
    #(Tint, lint, nodes_in_order, nodes_to_index, P, Q) = P1.make_unifrac_input_and_normalize(P2)
    #(unweighted, _) = EMDU.EMDUnifrac_unweighted(Tint, lint, nodes_in_order, P, Q)
    ###################################################################
    # compute normalizing factor attempt 1:
    # max unweighted unifrac is much more difficult: will need to take the ground truth \union other profile tree,
    # stick non-zero weights on all the nodes that don't have abundance on the ground truth part, then re-run it
    # through emdUnifrac.
    #Q_worst = np.ones(Q.size)
    #Q_worst[np.nonzero(P)[0]] = 0
    #(max_unweighted, _) = EMDU.EMDUnifrac_unweighted(Tint, lint, nodes_in_order, P, Q_worst)
    # note: this is the worst possible unweighted UniFrac on the gs \union other_profile tree. This actually
    #  leads to misleading results if that's what we normalize by, since a) technically the unweighted unifrac can
    #  be unbounded, and b) the results cannon be compared across different tools (since they all have different
    #  gs \union other_profile trees
    ###################################################################
    # compute normalizing factor attempt 2:
    # instead, compute the unweighted unifrac for the gold standard to an empty profile and the same for the
    # input profile, then take the (log?) ratio of these two

    ##(gs_unweighted_unifrac, _) = EMDU.EMDUnifrac_unweighted(Tint, lint, nodes_in_order, P, np.zeros(P.shape))
    ##(other_unweighted_unifrac, _) = EMDU.EMDUnifrac_unweighted(Tint, lint, nodes_in_order, Q, np.zeros(Q.shape))

    # just the ratio
    # unweighted = other_unweighted_unifrac / float(gs_unweighted_unifrac)
    # 1 - ratio (so zero is best and smaller is better)
    #unweighted = 1- other_unweighted_unifrac / float(gs_unweighted_unifrac)
    # relative difference

    ##unweighted = np.abs(other_unweighted_unifrac - gs_unweighted_unifrac) / float(gs_unweighted_unifrac)

    # Just take the log
    #unweighted = np.log(unweighted) if unweighted > 1 else 0
    ###################################################################
    # compute normalizing factor attempt 3:
    # ||WP - WQ|| / ||WP||
    # Basically: relative unweighted unifrac error. This is the best so far
    (gs_unweighted_unifrac, _) = EMDU.EMDUnifrac_unweighted(Tint, lint, nodes_in_order, P, np.zeros(P.shape))
    (raw_unweighted_unifrac, _) = EMDU.EMDUnifrac_unweighted(Tint, lint, nodes_in_order, P, Q)
    unweighted = raw_unweighted_unifrac / float(gs_unweighted_unifrac)

    return weighted, unweighted


def weighted_unifrac(pf1, pf2):
    P1 = copy.deepcopy(pf1)
    P2 = copy.deepcopy(pf2)
    (Tint, lint, nodes_in_order, nodes_to_index, P, Q) = P1.make_unifrac_input_and_normalize(P2)
    (val, _) = EMDU.EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
    return val


def unweighted_unifrac(pf1, pf2):
    P1 = copy.deepcopy(pf1)
    P2 = copy.deepcopy(pf2)
    (Tint, lint, nodes_in_order, nodes_to_index, P, Q) = P1.make_unifrac_input_and_normalize(P2)
    (val, _) = EMDU.EMDUnifrac_unweighted(Tint, lint, nodes_in_order, P, Q)
    return val


def print_list_unifrac(unifrac_list, labels, stream=sys.stdout):
    labels_iterator = iter(labels)
    stream.write("tool\tweighted\tunweighted\n")
    for unifrac in unifrac_list:
        stream.write("{}\t".format(next(labels_iterator)))
        stream.write("\t".join([format(unifrac[0], '.3f'), format(unifrac[1], '.3f')]))
        stream.write("\n")


def read_params(args):
    parser = argparse.ArgumentParser(description='')
    arg = parser.add_argument
    arg('--input', metavar='files_file', type=str, required=True,
        default=None, help="File of CAMI profile files to use")
    arg('--output', metavar='output_file', required=True, default=None, type=str,
        help="Output file (you should have this end in .csv as it is a matrix)")
    arg('--threshold', metavar='threshold', required=True, default=None, type=float,
        help="Value to threshold profiles to before computing EMDUnifrac. "
             "NOTE THIS VALUE IS IN PERCENTAGES so if you want 1% use 1")
    return vars(parser.parse_args())


if __name__ == '__main__':
    par = read_params(sys.argv)
    files_file = par['input']
    output_file = par['output']
    threshold = par['threshold']

    # Get all the profile file names
    files = []
    fid = open(files_file, 'r')
    for line in fid.readlines():
        files.append(line.strip())
    fid.close()

    # Import all the profiles
    profiles = []
    for file_name in files:
        profiles.append(PF.Profile(input_file_name=file_name))

    # Threshold all the profiles
    for profile in profiles:
        profile.threshold(threshold=threshold)

    # Compute EMDUnifrac
    D = np.zeros((len(profiles), len(profiles)))
    for i in range(len(profiles)):
        for j in range(i + 1, len(profiles)):
            (wu, uu) = compute_unifrac(profiles[i], profiles[j])
            D[i, j] = wu
            D[j, i] = uu

    # Save results in tsv
    np.savetxt(output_file, D, delimiter='\t', newline='\n')
