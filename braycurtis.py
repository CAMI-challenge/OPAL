from scipy.spatial import distance


def braycurtis(gs_rank_to_taxid_to_percentage, rank_to_taxid_to_percentage):
    rank_to_braycurtis = {}
    for rank in gs_rank_to_taxid_to_percentage:
        if rank not in rank_to_taxid_to_percentage:
            continue
        gs = []
        predictions = []
        # TODO: change to iterate over union
        # iterate over intersection of taxa
        for taxid in set(gs_rank_to_taxid_to_percentage[rank].keys()) & set(rank_to_taxid_to_percentage[rank].keys()):
            gs.append(gs_rank_to_taxid_to_percentage[rank][taxid])
            predictions.append(rank_to_taxid_to_percentage[rank][taxid])
        if len(gs) > 0:
            rank_to_braycurtis[rank] = distance.braycurtis(gs, predictions)
    return rank_to_braycurtis
