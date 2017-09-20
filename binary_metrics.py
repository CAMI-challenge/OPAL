#!/usr/bin/env python

"""
This script computes the following binary metrics
  * tp
  * fp
  * fn
  * precision
  * recall
  * Jaccard index

Tests can be executed by running
> python binary_metrics.py
"""


def __get_existing_taxa(rank):
    """ Return set of taxids with abundance > 0

    >>> __get_existing_taxa(query_rank)
    [123]

    :param rank: Set of taxids of specific rank
    :return: list of taxids
    """
    return list(k for k, v in rank.items() if v > 0)


def __get_non_existing_taxa(rank_query, rank_truth):
    """Return set of taxids with abundance <= 0
    >>> __get_non_existing_taxa(query_rank, truth_rank)
    []
    >>> __get_non_existing_taxa(query_rank2, truth_rank)
    [122]

    :param rank: Set of taxids of specific rank
    :return: list of taxids
    """
    rank_truth_taxids = set(k for k, v in rank_truth.items() if v > 0)
    rank_query_taxids = set(k for k, v in rank_query.items() if v > 0)
    return list(rank_query_taxids.difference(rank_truth_taxids))


def __get_tp(rank_query, rank_truth):
    """ Returns true positive
    >>> __get_tp(test_query_rank, test_truth_rank)
    1.0

    """
    rank_query_taxids = __get_existing_taxa(rank_query)
    rank_truth_taxids = __get_existing_taxa(rank_truth)
    return float(len(list(set(rank_query_taxids).intersection(rank_truth_taxids))))


def __get_fp(rank_query, rank_truth):
    """ Returns false positive
    >>> __get_fp(test_query_rank, test_truth_rank)
    0.0

    """

    rank_query_taxids = __get_existing_taxa(rank_query)
    rank_truth_taxids = __get_non_existing_taxa(rank_query, rank_truth)
    return float(len(list(set(rank_query_taxids).intersection(rank_truth_taxids))))


def __get_fn(rank_query, rank_truth):
    """ Returns false negative
    >>> __get_fn(test_query_rank, test_truth_rank)
    0.0

    """
    rank_query_taxids = __get_non_existing_taxa(rank_truth, rank_query)
    rank_truth_taxids = __get_existing_taxa(rank_truth)
    return float(len(list(set(rank_query_taxids).intersection(rank_truth_taxids))))


def precision(tp, fp):
    return tp / (tp + fp)


def recall(tp, fn):
    return tp / (tp + fn)


def jaccard_index(rank_query, rank_truth):
    """ Returns the Jaccard index
    >>> jaccard_index(test_query_rank, test_truth_rank)
    1.0

    """
    rank_query = __get_existing_taxa(rank_query)
    rank_truth = __get_existing_taxa(rank_truth)
    intersection = len(list(set(rank_query).intersection(rank_truth)))
    union = len(list(set(rank_query).union(rank_truth)))
    return intersection / union


def compute_rank_metrics(rank_query, rank_truth):
    """ Returns metrics for one rank
    >>> compute_rank_metrics(test_query_rank, test_truth_rank)
    (1.0, 1.0, 1.0, 0.0, 0.0, 1.0)

    """
    tp = __get_tp(rank_query, rank_truth)
    fn = __get_fn(rank_query, rank_truth)
    fp = __get_fp(rank_query, rank_truth)
    return precision(tp, fp), recall(tp, fn), tp, fn, fp, jaccard_index(rank_query, rank_truth)


def compute_tree_metrics(query, truth):
    """ Return metrics for tree
    >>> compute_tree_metrics(query_tree, truth_tree)
    {'species': (0.4, 0.6666666666666666, 2.0, 1.0, 3.0, 0.3333333333333333)}
    """
    def check_for_rank(query, rank):
        """Make sure that empty level is produced
        """
        if rank in query:
            return query[rank]
        else:
            return {}
    return {rank: compute_rank_metrics(check_for_rank(query, rank), taxids) for rank, taxids in truth.items()}


if __name__ == "__main__":
    import doctest

    test_query_rank = dict()
    test_query_rank[123] = 0.1
    test_query_rank[1232] = 0.0

    test_truth_rank = dict()
    test_truth_rank[123] = 0.1
    test_truth_rank[1232] = 0.0

    test_query_rank2 = dict()
    test_query_rank2[123] = 0.1
    test_query_rank2[122] = 4.0

    test_truth_tree = dict()
    test_truth_tree["species"] = { "A" : 0.5,
                                   "B" : 0.3,
                                   "E" : 0.7}

    test_query_tree = dict()
    test_query_tree["species"] = { "A" : 0.5,
                                   "B" : 0.9,
                                   "C" : 0.4,
                                   "D" : 0.2,
                                   "F" : 0.2}




    doctest.testmod(extraglobs={'query_rank': test_query_rank,
                                'query_rank2': test_query_rank2,
                                'truth_rank': test_truth_rank,
                                'truth_tree': test_truth_tree,
                                'query_tree': test_query_tree})
