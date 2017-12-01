from scipy.spatial import distance


def braycurtis(gs_rank_to_taxid_to_percentage, rank_to_taxid_to_percentage):
    rank_to_braycurtis = {}
    for rank in gs_rank_to_taxid_to_percentage:
        if rank not in rank_to_taxid_to_percentage:
            continue
        gs = []
        predictions = []
        taxa_union = set(gs_rank_to_taxid_to_percentage[rank].keys()).union(rank_to_taxid_to_percentage[rank].keys())
        for taxid in taxa_union:
            gs.append(gs_rank_to_taxid_to_percentage[rank][taxid] if taxid in gs_rank_to_taxid_to_percentage[rank] else .0)
            predictions.append(rank_to_taxid_to_percentage[rank][taxid] if taxid in rank_to_taxid_to_percentage[rank] else .0)
        if len(gs) > 0:
            rank_to_braycurtis[rank] = distance.braycurtis(gs, predictions)
    return rank_to_braycurtis
