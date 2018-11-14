#!/usr/bin/env python

import docker
import time
import argparse
import os
import errno


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
    parser = argparse.ArgumentParser()
    parser.add_argument('--image', help='Docker image of profiler', required=True)
    parser.add_argument('--yaml', help='Bioboxes YAML file', required=False)
    parser.add_argument('--input_dir', help='Input directory containing gzipped FASTQ files', required=True)
    parser.add_argument('--output_dir', help='Output directory', required=True)

    args = parser.parse_args()

    volumes = dict()
    volumes[args.input_dir] = {'bind': '/bbx/mnt/input', 'mode': 'ro'}
    yaml = None
    if args.yaml:
        if os.path.dirname(args.yaml) == args.input_dir:
            yaml = {'YAML': '/bbx/mnt/input/' + os.path.basename(args.yaml)}
        else:
            yaml = {'YAML': '/bbx/mnt/yaml/' + os.path.basename(args.yaml)}
            volumes[os.path.dirname(args.yaml)] = {'bind': '/bbx/mnt/yaml', 'mode': 'ro'}
    volumes[args.output_dir] = {'bind': '/bbx/mnt/output', 'mode': 'rw'}
    volumes[os.path.join(args.output_dir, 'metadata')] = {'bind': '/bbx/metadata', 'mode': 'rw'}
    volumes[os.path.join(args.output_dir, 'cache')] = {'bind': '/cache', 'mode': 'rw'}

    make_sure_path_exists(args.output_dir)
    make_sure_path_exists(os.path.join(args.output_dir, 'metadata'))
    make_sure_path_exists(os.path.join(args.output_dir, 'cache'))

    start_time = time.time()
    max_total_rss = run_docker(args.image, volumes, yaml)
    elapsed_time = time.time() - start_time
    memory = max_total_rss / 1048576.0

    print('{} MB\n{:.2f} seconds\n'.format(memory, elapsed_time))

    with open(os.path.join(args.output_dir, 'maxmemory_runtime.txt'), 'w') as f:
        f.write('{} MB\n{:.2f} seconds\n'.format(memory, elapsed_time))


if __name__ == "__main__":
    main()
