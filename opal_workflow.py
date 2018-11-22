#!/usr/bin/env python3

import argparse
import os
import re
import sys
import subprocess
import logging
from opal import get_labels
from opal import make_sure_path_exists
from opal_stats import get_image_dir_name
from version import __version__


def get_sorted_list_of_profiles_files(path):
    p = re.compile('^result_[0-9]+__.+\.profile$')
    all_files = os.listdir(path)
    mylist = []
    for file in all_files:
        if not p.match(file):
            continue
        mylist.append(file)
    mylist.sort()
    return mylist


def read_stats(results_dir, images_list):
    time = ''
    memory = ''
    comma = ''
    for image in images_list:
        image_dir_name = get_image_dir_name(image)
        with open(os.path.join(results_dir, image_dir_name, 'runtime_maxmemory.txt')) as f_input:
            # assumes first line is time and second, memory
            time += comma + str(float(f_input.readline().split(' ', 1)[0])/3600)
            memory += comma + str(float(f_input.readline().split(' ', 1)[0])/1024)
        comma = ','
    return time, memory


def preprocess_results(results_dir, images_list):
    index_empty_results = []
    p1 = re.compile('[0-9]+')
    p2 = re.compile('^@[^@]')
    for i, image in enumerate(images_list):
        image_dir_name = get_image_dir_name(image)
        file_path = os.path.join(results_dir, image_dir_name, 'all_results.profile')
        f_output = open(file_path, 'w')
        sorted_files = get_sorted_list_of_profiles_files(os.path.join(results_dir, image_dir_name))

        for file in sorted_files:
            sample_number = p1.search(file).group()
            f_input = open(os.path.join(results_dir, image_dir_name, file))

            for line in f_input:
                if p2.match(line):
                    key, value = line[1:].split(':', 1)
                    if key.upper() == 'SAMPLEID':
                        f_output.write('@{}:{}\n'.format(key, sample_number))
                    else:
                        f_output.write(line)
                else:
                    f_output.write(line)
            f_input.close()
        f_output.close()
        if os.stat(file_path).st_size == 0:
            index_empty_results.append(i)
    return index_empty_results


def get_profiles_list(results_dir, images_list):
    profiles_list = []
    for image in images_list:
        image_dir_name = get_image_dir_name(image)
        profiles_list.append(os.path.join(results_dir, image_dir_name, 'all_results.profile'))
    return profiles_list


def get_logger(output_dir):
    make_sure_path_exists(output_dir)
    logger = logging.getLogger('opal_workflow')
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    logging_fh = logging.FileHandler(os.path.join(output_dir, 'log.txt'))
    logging_fh.setFormatter(formatter)
    logger.addHandler(logging_fh)

    logging_stdout = logging.StreamHandler(sys.stdout)
    logging_stdout.setFormatter(formatter)
    logger.addHandler(logging_stdout)
    return logger


def main():
    parser = argparse.ArgumentParser(description='Run bioboxes of profilers with opal_stats and assess results with OPAL', add_help=False)
    group1 = parser.add_argument_group('required arguments')
    group1.add_argument('images', nargs='+', help='Docker images (bioboxes) of profilers')
    group1.add_argument('--input_dir', help='Input directory containing gzipped FASTQ files', required=True)
    group1.add_argument('--output_dir', help='Output directory', required=True)
    group1.add_argument('-g', '--gold_standard_file', help='Gold standard file', required=True)
    group2 = parser.add_argument_group('optional arguments')
    group2.add_argument('--yaml', help='Bioboxes YAML file (default: INPUT_DIR/biobox.yaml)', required=False)
    group2.add_argument("--volume", help='Docker volume', action='append')
    group2.add_argument('-n', '--no_normalization', help='Do not normalize samples', action='store_true')
    group2.add_argument('-p', '--plot_abundances', help='Plot abundances in the gold standard', action='store_true')
    group2.add_argument('-l', '--labels', help='Comma-separated names of profilers to be shown in OPAL', required=False)
    group2.add_argument('-d', '--desc', help='Description for HTML page', required=False)
    group2.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    group2.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    args = parser.parse_args()

    logger = get_logger(args.output_dir)
    images = args.images[:]
    labels = get_labels(args.labels, images)

    for image in images:
        logger.info('Running {}...'.format(image))
        parameters = [image,
                      '--input_dir=' + args.input_dir,
                      '--output_dir=' + args.output_dir]
        if args.yaml:
            parameters.append('--yaml=' + args.yaml)
        if args.volume:
            for volume in args.volume:
                parameters.append('--volume=' + volume)
        try:
            subprocess.run([sys.executable, 'opal_stats.py'] + parameters, check=True)
        except subprocess.CalledProcessError:
            logger.error('Error: opal_stats.py returned non-zero exit status for ' + image)
    logger.info('done')

    logger.info('Reading results...')
    index_empty_results = preprocess_results(args.output_dir, images)
    for index in sorted(index_empty_results, reverse=True):
        logger.warning('No results available for ' + images[index])
        del images[index]
        del labels[index]
    profiles_list = get_profiles_list(args.output_dir, images)
    time, memory = read_stats(args.output_dir, images)
    logger.info('done')

    logger.info('Running OPAL...')
    parameters = ['--gold_standard_file=' + args.gold_standard_file,
                  '--time=' + time,
                  '--memory=' + memory,
                  '--output_dir=' + os.path.join(args.output_dir, 'opal_output')]
    if args.no_normalization:
        parameters.append('--no_normalization')
    if args.plot_abundances:
        parameters.append('--plot_abundances')
    parameters.append('--labels=' + ','.join(labels))
    if args.desc:
        parameters.append('--desc=' + args.desc)
    parameters += profiles_list
    subprocess.run([sys.executable, 'opal.py'] + parameters, check=True)


if __name__ == "__main__":
    main()
