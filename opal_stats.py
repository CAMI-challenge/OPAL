#!/usr/bin/env python3

import docker
import time
import argparse
import os
import errno
from version import __version__


def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def run_docker(image, volumes, yaml):
    client = docker.from_env()
    container = client.containers.run(image=image, command=None, volumes=volumes, detach=True, environment=yaml)
    stats = container.stats(stream=True, decode=True)
    max_total_rss = .0

    while not container.status.startswith("e"):
        container.reload()
        next_stats = next(stats)
        time.sleep(1)
        try:
            total_rss = next_stats["memory_stats"]["stats"]["total_rss"]
            if total_rss > max_total_rss:
                max_total_rss = total_rss
        except KeyError:
            pass
    return max_total_rss


def main():
    parser = argparse.ArgumentParser(add_help=False)
    group1 = parser.add_argument_group('required arguments')
    group1.add_argument('image', nargs=1, help='Docker image (biobox) of profiler')
    group1.add_argument('--input_dir', help='Input directory containing gzipped FASTQ files', required=True)
    group1.add_argument('--output_dir', help='Output directory', required=True)
    group2 = parser.add_argument_group('optional arguments')
    group2.add_argument('--yaml', help='Bioboxes YAML file (default: INPUT_DIR/biobox.yaml)', required=False)
    group2.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    group2.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    args = parser.parse_args()
    input_dir = args.input_dir
    output_dir = args.output_dir
    image_dir_name = args.image.replace('/', '-')

    volumes = dict()
    yaml = None
    if args.yaml:
        if os.path.dirname(args.yaml) == input_dir:
            yaml = {'YAML': '/bbx/mnt/input/' + os.path.basename(args.yaml)}
        else:
            yaml = {'YAML': '/bbx/mnt/yaml/' + os.path.basename(args.yaml)}
            volumes[os.path.dirname(args.yaml)] = {'bind': '/bbx/mnt/yaml', 'mode': 'ro'}

    volumes[input_dir] = {'bind': '/bbx/mnt/input', 'mode': 'ro'}
    volumes[os.path.join(output_dir, image_dir_name)] = {'bind': '/bbx/mnt/output', 'mode': 'rw'}
    volumes[os.path.join(output_dir, image_dir_name, 'metadata')] = {'bind': '/bbx/metadata', 'mode': 'rw'}
    volumes[os.path.join(output_dir, image_dir_name, 'cache')] = {'bind': '/cache', 'mode': 'rw'}

    make_sure_path_exists(output_dir)
    make_sure_path_exists(os.path.join(output_dir, image_dir_name))
    make_sure_path_exists(os.path.join(output_dir, image_dir_name, 'metadata'))
    make_sure_path_exists(os.path.join(output_dir, image_dir_name, 'cache'))

    start_time = time.time()
    max_total_rss = run_docker(args.image, volumes, yaml)
    elapsed_time = time.time() - start_time
    memory = max_total_rss / 1048576.0

    print('{:.2f} seconds\n{} MB'.format(elapsed_time, memory))

    with open(os.path.join(output_dir, image_dir_name, 'runtime_maxmemory.txt'), 'w') as f:
        f.write('{:.2f} seconds\n{} MB\n'.format(elapsed_time, memory))


if __name__ == "__main__":
    main()
