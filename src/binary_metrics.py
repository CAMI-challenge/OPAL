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


class RankMetrics:
    """Class for saving metrics per rank
    """

    def __init__(self, rank):
        self.__rank = rank

    @property
    def rank(self):
        return self.__rank

    @property
    def tp(self):
        return self.__tp

    @property
    def fn(self):
        return self.__fn

    @property
    def fp(self):
        return self.__fp

    @property
    def precision(self):
        return self.__precision

    @property
    def recall(self):
        return self.__recall

    @property
    def jaccard(self):
        return self.__jaccard

    @property
    def f1(self):
        return self.__f1

    @rank.setter
    def rank(self, rank):
        self.__rank = rank

    @fn.setter
    def fn(self, fn):
        self.__fn = fn

    @tp.setter
    def tp(self, tp):
        self.__tp = tp

    @fp.setter
    def fp(self, fp):
        self.__fp = fp

    @precision.setter
    def precision(self, precision):
        self.__precision = precision

    @recall.setter
    def recall(self, recall):
        self.__recall = recall

    @jaccard.setter
    def jaccard(self, jaccard):
        self.__jaccard = jaccard

    @f1.setter
    def f1(self, f1):
        self.__f1 = f1

    def get_ordered_dict(self):
        from collections import OrderedDict
        return OrderedDict(sorted(self.__dict__.items()))

    def get_dict(self):
        return self.__dict__

    def get_pretty_dict(self):
        return {metric.split("_")[3]: value for metric, value in self.__dict__.items()}


def __get_taxa(rank):
    """ Return set of taxids with abundance > 0

    >>> __get_taxa(query_rank)
    {123}

    :param rank: Set of taxids of specific rank
    :return: list of taxids
    """
    return set(k for k, v in rank.items() if v > 0)


def __get_tp(rank_query_taxids, rank_truth_taxids):
    """ Returns true positive
    >>> __get_tp(test_rank_query_taxids, test_rank_truth_taxids)
    2

    """
    return len(rank_query_taxids.intersection(rank_truth_taxids))


def __get_fp(rank_query_fp_taxids):
    """ Returns false positive
    >>> __get_fp(test_rank_query_fp_taxids)
    0

    """
    return len(rank_query_fp_taxids)


def __get_fn(rank_query_fn_taxids):
    """ Returns false negative
    >>> __get_fn(test_rank_query_fn_taxids)
    0

    """
    return len(rank_query_fn_taxids)


def precision(tp, fp):
    if (tp + fp) == 0:
        return 0
    return tp / (tp + fp)


def recall(tp, fn):
    if (tp + fn) == 0:
        return 0
    return tp / (tp + fn)


def jaccard_index(tp, rank_query_taxids, rank_truth_taxids):
    """ Returns the Jaccard index
    >>> jaccard_index(test_tp, test_rank_query_taxids, test_rank_truth_taxids)
    1.0

    """
    union = len(rank_query_taxids.union(rank_truth_taxids))
    if union > 0:
        return tp / union
    else:
        return .0


def f1_score(this_precision, this_recall):
    """ Returns f1 score
    >>> f1_score(test_rank_metrics.precision, test_rank_metrics.recall)
    1.0

    """
    if this_precision + this_recall > 0:
        return 2 * this_precision * this_recall / (this_precision + this_recall)
    else:
        return None


def compute_rank_metrics(rank_query, rank_truth, rank):
    """ Returns metrics for one rank
    >>> compute_rank_metrics(test_query_rank, test_truth_rank, "species").get_ordered_dict()
    OrderedDict([('_RankMetrics__f1', 1.0), ('_RankMetrics__fn', 0), ('_RankMetrics__fp', 0), ('_RankMetrics__jaccard', 1.0), ('_RankMetrics__precision', 1.0), ('_RankMetrics__rank', 'species'), ('_RankMetrics__recall', 1.0), ('_RankMetrics__tp', 1)])

    """
    rank_query_taxids = __get_taxa(rank_query)
    rank_truth_taxids = __get_taxa(rank_truth)
    rank_query_fp_taxids = rank_query_taxids.difference(rank_truth_taxids)
    rank_query_fn_taxids = rank_truth_taxids.difference(rank_query_taxids)

    rank_metrics = RankMetrics(rank)
    rank_metrics.tp = __get_tp(rank_query_taxids, rank_truth_taxids)
    rank_metrics.fn = __get_fn(rank_query_fn_taxids)
    rank_metrics.fp = __get_fp(rank_query_fp_taxids)

    rank_metrics.recall = recall(rank_metrics.tp, rank_metrics.fn)
    rank_metrics.precision = precision(rank_metrics.tp, rank_metrics.fp)
    rank_metrics.f1 = f1_score(rank_metrics.precision, rank_metrics.recall)
    rank_metrics.jaccard = jaccard_index(rank_metrics.tp, rank_query_taxids, rank_truth_taxids)

    return rank_metrics


def compute_tree_metrics(query, truth):
    """ Return metrics for tree
    >>> compute_tree_metrics(query_tree, truth_tree)["species"].get_ordered_dict()
    OrderedDict([('_RankMetrics__f1', 0.5), ('_RankMetrics__fn', 1), ('_RankMetrics__fp', 3), ('_RankMetrics__jaccard', 0.3333333333333333), ('_RankMetrics__precision', 0.4), ('_RankMetrics__rank', 'species'), ('_RankMetrics__recall', 0.6666666666666666), ('_RankMetrics__tp', 2)])
    """

    def check_for_rank(query, rank):
        """Make sure that empty level is produced
        """
        if rank in query:
            return query[rank]
        else:
            return {}

    return {rank: compute_rank_metrics(check_for_rank(query, rank), taxids, rank) for rank, taxids in truth.items()}


def print_all_metrics(tree_metrics, path):
    import csv

    def get_header(tree):
        not_empty_ranks = [metrics.get_pretty_dict() for rank, metrics in tree.items() if
                           bool(metrics.get_pretty_dict())]
        if len(not_empty_ranks) > 0:
            return not_empty_ranks[0].keys()
        else:
            return []

    with open(path, 'w') as metrcs_out:
        writer = csv.DictWriter(metrcs_out, fieldnames=get_header(tree_metrics),
                                restval='ignore', delimiter='\t')
        writer.writeheader()
        {rank: writer.writerow(metrics.get_pretty_dict()) for rank, metrics in tree_metrics.items()}


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

    test_rank_query_taxids = set(test_query_rank.keys())
    test_rank_truth_taxids = set(test_truth_rank.keys())
    test_rank_query_fp_taxids = test_rank_query_taxids.difference(test_rank_truth_taxids)
    test_rank_query_fn_taxids = test_rank_truth_taxids.difference(test_rank_query_taxids)
    test_tp = 2

    test_truth_tree = dict()
    test_truth_tree["species"] = {"A": 0.5,
                                  "B": 0.3,
                                  "E": 0.7}

    test_query_tree = dict()
    test_query_tree["species"] = {"A": 0.5,
                                  "B": 0.9,
                                  "C": 0.4,
                                  "D": 0.2,
                                  "F": 0.2}

    tree_metrics = dict()
    test_rank_metrics = RankMetrics("species")
    test_rank_metrics.tp = 0
    test_rank_metrics.fn = 0
    test_rank_metrics.fp = 1
    test_rank_metrics.precision = 1
    test_rank_metrics.recall = 1
    test_rank_metrics.jaccard = 1
    test_rank_metrics.f1 = 1
    tree_metrics["species"] = test_rank_metrics

    doctest.testmod(extraglobs={'query_rank': test_query_rank,
                                'query_rank2': test_query_rank2,
                                'truth_rank': test_truth_rank,
                                'truth_tree': test_truth_tree,
                                'tree_metrics': tree_metrics,
                                'query_tree': test_query_tree})
