#!/usr/bin/env python

import math
from collections import OrderedDict
from src.utils import constants as c


class Shannon:
    def __init__(self, rank):
        self.__rank = rank

    @property
    def rank(self):
        return self.__rank

    @property
    def diversity(self):
        return self.__diversity

    @property
    def equitability(self):
        return self.__equitability

    @rank.setter
    def rank(self, rank):
        self.__rank = rank

    @diversity.setter
    def diversity(self, diversity):
        self.__diversity = diversity

    @equitability.setter
    def equitability(self, equitability):
        self.__equitability = equitability

    def get_dict(self):
        return self.__dict__

    def get_pretty_dict(self):
        return {metric.split("_")[3]: value for metric, value in self.__dict__.items()}


# Compute the Shannon index (diversity and equitability) with log base e
# The higher the index, the more equally distributed the taxa are
def compute_shannon_index(rank_to_taxid_to_percentage):
    rank_to_shannon = OrderedDict()
    for rank in c.ALL_RANKS:
        if rank not in rank_to_taxid_to_percentage:
            continue
        rank_to_shannon[rank] = Shannon(rank)
        rank_to_shannon[rank].diversity = .0
        rank_to_shannon[rank].equitability = .0
        for taxid in rank_to_taxid_to_percentage[rank]:
            percentage = rank_to_taxid_to_percentage[rank][taxid] / 100.0
            if percentage > .0:
                rank_to_shannon[rank].diversity -= percentage * math.log(percentage)
        if len(rank_to_taxid_to_percentage[rank]) > 1:
            rank_to_shannon[rank].equitability = rank_to_shannon[rank].diversity / math.log(len(rank_to_taxid_to_percentage[rank]))
    return rank_to_shannon
