#!/usr/bin/env python

import docker
import time
import argparse


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
    parser.add_argument('-i', '--image', help='Docker image of profiler', required=True)
    parser.add_argument('-v', '--volume', type=lambda x: x.split(':'), action='append', required=True)
    parser.add_argument('-y', '--yaml', help='Bioboxes YAML file', required=True)
    args = parser.parse_args()

    if args.output_file and not args.label:
        print("Option --output_file requires --label")
        exit(2)

    volumes = dict()
    for i in range(0, len(args.volume)):
        if len(args.volume[0]) < 3:
            print("Missing required volumes")
            exit(2)
        volumes[args.volume[i][0]] = {'bind': args.volume[i][1], 'mode': args.volume[i][2]}

    yaml = None
    if args.yaml:
        yaml = {'YAML': args.yaml}

    start_time = time.time()
    max_total_rss = run_docker(args.image, volumes, yaml)
    elapsed_time = time.time() - start_time
    memory = max_total_rss / 1048576.0
    print("{} MB".format(memory))
    print("{:.2f} seconds".format(elapsed_time))


if __name__ == "__main__":
    main()
