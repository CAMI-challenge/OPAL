#!/usr/bin/env python

import biom
import sys
import os


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


def get_column_indices(column_name_to_index):
    if "TAXID" not in column_name_to_index:
        raise RuntimeError("Column not found: {}".format("TAXID"))
    if "RANK" not in column_name_to_index:
        raise RuntimeError("Column not found: {}".format("RANK"))
    if "PERCENTAGE" not in column_name_to_index:
        raise RuntimeError("Column not found: {}".format("PERCENTAGE"))
    if "TAXPATH" not in column_name_to_index:
        raise RuntimeError("Column not found: {}".format("TAXPATH"))
    index_taxid = column_name_to_index["TAXID"]
    index_rank = column_name_to_index["RANK"]
    index_percentage = column_name_to_index["PERCENTAGE"]
    index_taxpath = column_name_to_index["TAXPATH"]
    if "TAXPATHSN" in column_name_to_index:
        index_taxpathsn = column_name_to_index["TAXPATHSN"]
    else:
        index_taxpathsn = None
    return index_rank, index_taxid, index_percentage, index_taxpath, index_taxpathsn


def open_profile_from_tsv(file_path):
    header = {}
    new_header = {}
    column_name_to_index = {}
    profile = []
    samples_list = []

    with open(file_path) as read_handler:
        for line in read_handler:
            if len(line.strip()) == 0 or line.startswith("#"):
                continue
            line = line.rstrip('\n')

            if line.startswith("@@"):
                header = new_header
                new_header = {}

                for index, column_name in enumerate(line[2:].split('\t')):
                    column_name_to_index[column_name] = index

                # store profile for sample and empty profile array
                if 'SAMPLEID' in header and 'VERSION' in header and 'RANKS' in header:
                    if len(profile) > 0:
                        samples_list.append((header['SAMPLEID'], header, profile))
                        profile = []
                else:
                    print("Header in file {} is incomplete. Check if the header of each sample contains at least SAMPLEID, VERSION, and RANKS.".format(file_path))
                    raise
                continue

            # parse header
            if line.startswith("@"):
                key, value = line[1:].split(':', 1)
                new_header[key.upper()] = value.strip()
                continue

            index_rank, index_taxid, index_percentage, index_taxpath, index_taxpathsn = get_column_indices(column_name_to_index)
            row_data = line.split('\t')

            prediction = Prediction()
            prediction.taxid = row_data[index_taxid]
            prediction.rank = row_data[index_rank]
            prediction.percentage = float(row_data[index_percentage])
            prediction.taxpath = row_data[index_taxpath]
            if isinstance(index_taxpathsn, int):
                prediction.taxpathsn = row_data[index_taxpathsn]
            else:
                prediction.taxpathsn = None
            profile.append(prediction)

    # store profile for last sample
    if 'SAMPLEID' in header and 'VERSION' in header and 'RANKS' in header:
        if len(profile) > 0:
            samples_list.append((header['SAMPLEID'], header, profile))
    else:
        print("Header in file {} is incomplete. Check if the header of each sample contains at least SAMPLEID, VERSION, and RANKS.".format(file_path))
        raise

    return samples_list


def open_profile(file_path):
    if not os.path.exists(file_path):
        sys.exit("Input file {} does not exist.".format(file_path))
    try:
        table = biom.load_table(file_path)
    except:
        try:
            return open_profile_from_tsv(file_path)
        except:
            sys.exit("Input file could not be read.")

    samples_list = []
    samples = table.ids(axis='sample')
    ids = table.ids(axis='observation')

    for sample_id in samples:
        sample_metadata = table.metadata(id=sample_id, axis='sample')
        profile = []
        for id in ids:
            prediction = Prediction()
            metadata = table.metadata(id=id, axis='observation')
            prediction.taxid = id
            prediction.rank = metadata['rank']
            prediction.percentage = float(table.get_value_by_ids(obs_id=id, samp_id=sample_id))
            prediction.taxpath = metadata['taxpath']
            prediction.taxpathsn = metadata['taxpathsn']
            profile.append(prediction)
        samples_list.append((sample_id, sample_metadata, profile))
    return samples_list


def get_rank_to_taxid_to_percentage(profile):
    rank_to_taxid_to_percentage = {}
    for prediction in profile:
        if prediction.rank not in rank_to_taxid_to_percentage:
            rank_to_taxid_to_percentage[prediction.rank] = {}
        rank_to_taxid_to_percentage[prediction.rank][prediction.taxid] = prediction.percentage
    return rank_to_taxid_to_percentage
