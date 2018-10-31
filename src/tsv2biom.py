#!/usr/bin/env python

import sys
import argparse
import biom
from biom.util import biom_open
from collections import OrderedDict
from src.utils import load_data


def convert_to_biom(file_paths, output_file, is_json=False):
    tax_ids_row_mapping = OrderedDict()
    sample_ids = []
    sample_metadata = []
    data = {}
    observ_metadata = {}
    ids = []

    for file_path in file_paths:
        ids.append(file_path.split('/')[-1])
        try:
            samples_list = load_data.open_profile_from_tsv(file_path, False)
        except:
            sys.exit("Input file could not be read.")

        for i, sample in enumerate(samples_list):
            sample_id, header, profile = sample

            if sample_id == '':
                sample_id = 'S{}'.format(i)
                header['SAMPLEID'] = sample_id

            sample_ids.append(sample_id)
            sample_metadata.append(header)

            for prediction in profile:
                if prediction.taxid not in tax_ids_row_mapping:
                    tax_ids_row_mapping[prediction.taxid] = len(tax_ids_row_mapping)
                data[(tax_ids_row_mapping[prediction.taxid], i)] = prediction.percentage
                observ_metadata[prediction.taxid] = prediction.get_metadata()

    # sort metadata
    observ_metadata = [observ_metadata[key] for key in tax_ids_row_mapping.keys()]

    table = biom.Table(data=data,
                       observation_ids=list(tax_ids_row_mapping.keys()),
                       sample_ids=sample_ids,
                       observation_metadata=observ_metadata,
                       sample_metadata=sample_metadata,
                       table_id='_'.join(ids),
                       type='OTU table')

    if is_json:
        with open(output_file, 'w') as f:
            f.write(table.to_json(generated_by='OPAL'))
    else:
        with biom_open(output_file, 'w') as f:
            table.to_hdf5(f, 'OPAL')


def main():
    parser = argparse.ArgumentParser(description="Convert profile in CAMI Bioboxes format to BIOM format")
    parser.add_argument('files', nargs='+', help="Input file(s), one file per sample")
    parser.add_argument('-o', '--output_file', help="Output file", required=True)
    parser.add_argument('-j', '--json', action='store_true', help="Output in json (default: hdf5)")
    args = parser.parse_args()
    convert_to_biom(args.files, args.output_file, args.json)


if __name__ == "__main__":
    main()