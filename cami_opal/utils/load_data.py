#!/usr/bin/env python

import os
import logging
from collections import defaultdict


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
        logging.getLogger('opal').critical("Column not found: {}".format("TAXID"))
        raise RuntimeError
    if "RANK" not in column_name_to_index:
        logging.getLogger('opal').critical("Column not found: {}".format("RANK"))
        raise RuntimeError
    if "PERCENTAGE" not in column_name_to_index:
        logging.getLogger('opal').critical("Column not found: {}".format("PERCENTAGE"))
        raise RuntimeError
    if "TAXPATH" not in column_name_to_index:
        logging.getLogger('opal').critical("Column not found: {}".format("TAXPATH"))
        raise RuntimeError
    index_taxid = column_name_to_index["TAXID"]
    index_rank = column_name_to_index["RANK"]
    index_percentage = column_name_to_index["PERCENTAGE"]
    index_taxpath = column_name_to_index["TAXPATH"]
    if "TAXPATHSN" in column_name_to_index:
        index_taxpathsn = column_name_to_index["TAXPATHSN"]
    else:
        index_taxpathsn = None
    return index_rank, index_taxid, index_percentage, index_taxpath, index_taxpathsn


def normalize_samples(samples_list):
    for sample in samples_list:
        sample_id, sample_metadata, profile = sample
        sum_per_rank = defaultdict(float)
        for prediction in profile:
            sum_per_rank[prediction.rank] += prediction.percentage
        for prediction in profile:
            if prediction.percentage > 0:
                prediction.percentage = (prediction.percentage / sum_per_rank[prediction.rank]) * 100.0


def open_profile_from_tsv(file_path, normalize):
    header = {}
    column_name_to_index = {}
    profile = []
    samples_list = []
    predictions_dict = {}
    reading_data = False
    got_column_indices = False

    with open(file_path) as read_handler:
        for line in read_handler:
            if len(line.strip()) == 0 or line.startswith("#"):
                continue
            line = line.rstrip('\n')

            # parse header with column indices
            if line.startswith("@@"):
                for index, column_name in enumerate(line[2:].split('\t')):
                    column_name_to_index[column_name] = index
                index_rank, index_taxid, index_percentage, index_taxpath, index_taxpathsn = get_column_indices(column_name_to_index)
                got_column_indices = True
                reading_data = False
                continue

            # parse header with metadata
            if line.startswith("@"):
                # if last line contained sample data and new header starts, store profile for sample
                if reading_data:
                    if 'SAMPLEID' in header and 'VERSION' in header and 'RANKS' in header:
                        if len(profile) > 0:
                            samples_list.append((header['SAMPLEID'], header, profile))
                            profile = []
                            predictions_dict = {}
                    else:
                        logging.getLogger('opal').critical("Header in file {} is incomplete. Check if the header of each sample contains at least SAMPLEID, VERSION, and RANKS.\n".format(file_path))
                        raise RuntimeError
                    header = {}
                reading_data = False
                got_column_indices = False
                key, value = line[1:].split(':', 1)
                header[key.upper()] = value.strip()
                continue

            if not got_column_indices:
                logging.getLogger('opal').critical("Header line starting with @@ in file {} is missing or at wrong position.\n".format(file_path))
                raise RuntimeError

            reading_data = True
            row_data = line.split('\t')

            taxid = row_data[index_taxid]
            # if there is already a prediction for taxon, only sum abundance
            if taxid in predictions_dict:
                prediction = predictions_dict[taxid]
                prediction.percentage += float(row_data[index_percentage])
            else:
                if float(row_data[index_percentage]) == .0:
                    continue
                prediction = Prediction()
                predictions_dict[taxid] = prediction
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
        if reading_data and len(profile) > 0:
            samples_list.append((header['SAMPLEID'], header, profile))
    else:
        logging.getLogger('opal').critical("Header in file {} is incomplete. Check if the header of each sample contains at least SAMPLEID, VERSION, and RANKS.\n".format(file_path))
        raise RuntimeError

    if normalize:
        normalize_samples(samples_list)

    return samples_list


def open_profile(file_path, normalize):
    if not os.path.exists(file_path):
        logging.getLogger('opal').critical("Input file {} does not exist.".format(file_path))
        exit(1)
    try:
        return open_profile_from_tsv(file_path, normalize)
    except:
        logging.getLogger('opal').critical("Input file could not be read.")
        exit(1)

    # try:
    #     table = biom.load_table(file_path)
    # except:
    #     try:
    #         return open_profile_from_tsv(file_path, normalize)
    #     except:
    #         logging.getLogger('opal').critical("Input file could not be read.")
    #         exit(1)
    #
    # samples_list = []
    # samples = table.ids(axis='sample')
    # obs_ids = table.ids(axis='observation')
    #
    # for sample_id in samples:
    #     sample_metadata = table.metadata(id=sample_id, axis='sample')
    #     profile = []
    #     for obs_id in obs_ids:
    #         percentage = float(table.get_value_by_ids(obs_id=obs_id, samp_id=sample_id))
    #         if percentage == 0:
    #             continue
    #         prediction = Prediction()
    #         metadata = table.metadata(id=obs_id, axis='observation')
    #         taxpathsn_list = [x for x in metadata['taxonomy'] if x[3:]]
    #         taxpathsn = '|'.join(taxpathsn_list)
    #         prediction.taxid = taxpathsn_list[-1]
    #         prediction.rank = c.ALL_RANKS[len(taxpathsn_list) - 1]
    #         prediction.percentage = percentage
    #         prediction.taxpath = taxpathsn
    #         prediction.taxpathsn = taxpathsn
    #         profile.append(prediction)
    #     samples_list.append((sample_id, sample_metadata, profile))
    #
    # if normalize:
    #     normalize_samples(samples_list)
    #
    # return samples_list


def load_profiles(gold_standard_file, profiles_files, normalize):
    gs_samples_list = open_profile(gold_standard_file, normalize)
    sample_ids_list = []
    for sample in gs_samples_list:
        sample_id, sample_metadata, profile = sample
        sample_ids_list.append(sample_id)

    profiles_list_to_samples_list = []
    for profile_file in profiles_files:
        profiles_list_to_samples_list.append(open_profile(profile_file, normalize))

    return sample_ids_list, gs_samples_list, profiles_list_to_samples_list


def get_taxa_names(profile):
    tax_id_to_name = {}
    for prediction in profile:
        tax_id_to_name[prediction.taxid] = prediction.taxpathsn.split("|")[-1]
    return tax_id_to_name


def get_rank_to_taxid_to_percentage(profile, rank=None):
    rank_to_taxid_to_percentage = {}
    for prediction in profile:
        if rank and prediction.rank != rank:
            continue
        if prediction.rank not in rank_to_taxid_to_percentage:
            rank_to_taxid_to_percentage[prediction.rank] = {}
        rank_to_taxid_to_percentage[prediction.rank][prediction.taxid] = prediction.percentage
    return rank_to_taxid_to_percentage


def get_rank_to_taxid_to_percentage_filtered(rank_to_taxid_to_percentage, filter_tail_percentage):
    rank_to_taxid_to_percentage_filtered = defaultdict(lambda: defaultdict(float))
    for rank in rank_to_taxid_to_percentage:
        sorted_profile = sorted(rank_to_taxid_to_percentage[rank].items(), key=lambda x: x[1])
        sum_all = .0
        for item in sorted_profile:
            sum_all += item[1]
            if sum_all > filter_tail_percentage and item[1] > 0:
                rank_to_taxid_to_percentage_filtered[rank][item[0]] = item[1]
    return rank_to_taxid_to_percentage_filtered
