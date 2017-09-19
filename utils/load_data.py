#!/usr/bin/env python


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


def read_rows(input_stream, index_rank, index_taxid, index_percentage):
    for line in input_stream:
        if len(line.strip()) == 0 or line.startswith("#"):
            continue
        line = line.rstrip('\n')
        row_data = line.split('\t')
        rank = row_data[index_rank]
        taxid = row_data[index_taxid]
        percentage = row_data[index_percentage]
        yield rank, taxid, percentage


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
    return index_rank, index_taxid, index_percentage


def open_profile(file_path):
    rank_to_taxid_to_percentage = {}
    with open(file_path) as read_handler:
        index_rank, index_taxid, index_percentage = get_column_indices(read_handler)
        for rank, taxid, percentage in read_rows(read_handler, index_rank, index_taxid, index_percentage):
            if rank not in rank_to_taxid_to_percentage:
                rank_to_taxid_to_percentage[rank] = {}
            rank_to_taxid_to_percentage[rank][taxid] = float(percentage)
    return rank_to_taxid_to_percentage
