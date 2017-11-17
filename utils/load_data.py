#!/usr/bin/env python

import biom
import sys


class Prediction:
    def __init__(self):
        pass

    @property
    def rank(self):
        return self.__rank

    @property
    def taxid(self):
        return self.__taxid

    @property
    def percentage(self):
        return self.__percentage

    @property
    def taxpath(self):
        return self.__taxpath

    @property
    def taxpathsn(self):
        return self.__taxpathsn

    @rank.setter
    def rank(self, rank):
        self.__rank = rank

    @taxid.setter
    def taxid(self, taxid):
        self.__taxid = taxid

    @percentage.setter
    def percentage(self, percentage):
        self.__percentage = percentage

    @taxpath.setter
    def taxpath(self, taxpath):
        self.__taxpath = taxpath

    @taxpathsn.setter
    def taxpathsn(self, taxpathsn):
        self.__taxpathsn = taxpathsn

    def get_dict(self):
        return self.__dict__

    def get_pretty_dict(self):
        return {property.split("_")[3]: value for property, value in self.__dict__.items()}

    def get_metadata(self):
        return {'rank': self.__rank, 'taxpath': self.__taxpath, 'taxpathsn': self.__taxpathsn}


def read_header(input_stream):
    header = {}
    column_names = {}
    for line in input_stream:
        if len(line.strip()) == 0 or line.startswith("#"):
            continue
        line = line.rstrip('\n')
        if line.startswith("@@"):
            for index, column_name in enumerate(line[2:].split('\t')):
                column_names[column_name] = index
            return header, column_names
        if line.startswith("@"):
            key, value = line[1:].split(':', 1)
            header[key] = value.strip()


def read_rows(input_stream, index_rank, index_taxid, index_percentage, index_taxpath, index_taxpathsn):
    for line in input_stream:
        if len(line.strip()) == 0 or line.startswith("#"):
            continue
        line = line.rstrip('\n')
        row_data = line.split('\t')
        rank = row_data[index_rank]
        taxid = row_data[index_taxid]
        percentage = row_data[index_percentage]
        taxpath = row_data[index_taxpath]
        taxpathsn = row_data[index_taxpathsn]
        yield rank, taxid, percentage, taxpath, taxpathsn


def get_column_indices(input_stream):
    header, column_names = read_header(input_stream)
    if "TAXID" not in column_names:
        raise RuntimeError("Column not found: {}".format("TAXID"))
    if "RANK" not in column_names:
        raise RuntimeError("Column not found: {}".format("RANK"))
    if "PERCENTAGE" not in column_names:
        raise RuntimeError("Column not found: {}".format("PERCENTAGE"))
    index_taxid = column_names["TAXID"]
    index_rank = column_names["RANK"]
    index_percentage = column_names["PERCENTAGE"]
    if "TAXPATH" in column_names:
        index_taxpath = column_names["TAXPATH"]
    if "TAXPATHSN" in column_names:
        index_taxpathsn = column_names["TAXPATHSN"]
    return header, index_rank, index_taxid, index_percentage, index_taxpath, index_taxpathsn


def open_profile_from_tsv(file_path):
    profile = []
    with open(file_path) as read_handler:
        header, index_rank, index_taxid, index_percentage, index_taxpath, index_taxpathsn = get_column_indices(read_handler)
        for rank, taxid, percentage, taxpath, taxpathsn in read_rows(read_handler, index_rank, index_taxid, index_percentage, index_taxpath, index_taxpathsn):
            prediction = Prediction()
            prediction.taxid = taxid
            prediction.rank = rank
            prediction.percentage = float(percentage)
            prediction.taxpath = taxpath
            prediction.taxpathsn = taxpathsn
            profile.append(prediction)
    return header, profile


def open_profile(file_path):
    try:
        table = biom.load_table(file_path)
    except:
        try:
            return open_profile_from_tsv(file_path)
        except:
            sys.exit("Incorrect file format of input profile")

    samples = table.ids(axis='sample')
    ids = table.ids(axis='observation')

    for sample in samples:
        sample_metadata = table.metadata(id=sample, axis='sample')
        profile = []
        for id in ids:
            prediction = Prediction()
            metadata = table.metadata(id=id, axis='observation')
            prediction.taxid = id
            prediction.rank = metadata['rank']
            prediction.percentage = float(table.get_value_by_ids(obs_id=id, samp_id=sample))
            prediction.taxpath = metadata['taxpath']
            prediction.taxpathsn = metadata['taxpathsn']
            profile.append(prediction)
        # TODO: return profile of multiple samples
        return sample_metadata, profile


def get_rank_to_taxid_to_percentage(profile):
    rank_to_taxid_to_percentage = {}
    for prediction in profile:
        if prediction.rank not in rank_to_taxid_to_percentage:
            rank_to_taxid_to_percentage[prediction.rank] = {}
        rank_to_taxid_to_percentage[prediction.rank][prediction.taxid] = prediction.percentage
    return rank_to_taxid_to_percentage
