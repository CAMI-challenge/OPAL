#!/usr/bin/env python3

import argparse
import os
import re
import sys
import subprocess


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
        image_dir_name = image.replace('/', '-')
        with open(os.path.join(results_dir, image_dir_name, 'runtime_maxmemory.txt')) as f_input:
            # assumes first line is time and second, memory
            time += comma + f_input.readline().split(' ', 1)[0]
            memory += comma + f_input.readline().split(' ', 1)[0]
        comma = ','
    return time, memory


def preprocess_results(results_dir, images_list):
    p1 = re.compile('[0-9]+')
    p2 = re.compile('^@[^@]')
    for image in images_list:
        image_dir_name = image.replace('/', '-')
        f_output = open(os.path.join(results_dir, image_dir_name, 'all_results.profile'), 'w')
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


def get_profiles_list(results_dir, images_list):
    profiles_list = []
    for image in images_list:
        image_dir_name = image.replace('/', '-')
        profiles_list.append(os.path.join(results_dir, image_dir_name, 'all_results.profile'))
    return profiles_list


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gold_standard_file', help='Gold standard file', required=True)
    parser.add_argument('--images', help='Comma-separated list of Docker images of profilers', required=True)
    parser.add_argument('--yaml', help='Bioboxes YAML file', required=False)
    parser.add_argument('--input_dir', help='Input directory containing gzipped FASTQ files', required=True)
    parser.add_argument('--output_dir', help='Output directory', required=True)
    parser.add_argument('-d', '--desc', help='Description for HTML page (optional)', required=False)
    parser.add_argument('-l', '--labels', help='Comma-separated names of profilers to be shown in OPAL (optional)', required=False)
    args = parser.parse_args()

    images_list = [x.strip() for x in args.images.split(',')]

    for image in images_list:
        image_dir_name = image.replace('/', '-')
        parameters = ['--image=' + image,
                      '--yaml=' + args.yaml,
                      '--input_dir=' + args.input_dir,
                      '--output_dir=' + os.path.join(args.output_dir, image_dir_name)]
        try:
            subprocess.run([sys.executable, 'opal_stats.py'] + parameters, check=True)
        except subprocess.CalledProcessError:
            print('Error: opal_stats.py returned non-zero exit status 1 for ' + image, file=sys.stderr)

    preprocess_results(args.output_dir, images_list)
    profiles_list = get_profiles_list(args.output_dir, images_list)
    time, memory = read_stats(args.output_dir, images_list)

    parameters = ['--gold_standard_file=' + args.gold_standard_file,
                  '--time=' + time,
                  '--memory=' + memory,
                  '--output_dir=' + os.path.join(args.output_dir, 'opal_output')]
    if args.labels:
        parameters.append('--labels=' + args.labels)
    if args.desc:
        parameters.append('--desc=' + args.desc)
    parameters += profiles_list
    subprocess.run([sys.executable, 'opal.py'] + parameters, check=True)


if __name__ == "__main__":
    main()
