[![CircleCI](https://circleci.com/gh/CAMI-challenge/OPAL.svg?style=svg)](https://circleci.com/gh/CAMI-challenge/OPAL)

# OPAL
Profiling Assessment

# Requirements

* python &ge; 3.5
* numpy &ge; 1.13.0
* matplotlib &ge; 2.0.2
* dendropy &ge; 4.3.0
* pandas &ge; 0.20.3

# User Guide

## Installation

Download OPAL:
~~~BASH
wget https://github.com/CAMI-challenge/OPAL/archive/master.zip -O opal.zip
~~~
Or clone it using git:
~~~BASH
git clone https://github.com/CAMI-challenge/OPAL.git
~~~

Install dependencies as follows (tested on Linux Ubuntu 16.04):

~~~BASH
sudo apt-get install python3-pip
cd OPAL/
pip3 install -r requirements/default.txt --user
~~~

## Input
OPAL uses at least two files:
1. A gold standard taxonomic profile
2. One or more taxonomic profiles to be assessed

Files must be in the [CAMI profiling Bioboxes format](https://github.com/bioboxes/rfc/tree/master/data-format).

## Computed metrics

* Unifrac error
* L1 norm error
* True positives, false positives, false negatives
* Precision
* Recall
* F1 score
* Jaccard index
* Shannon diversity and equitability indices

## Running _opal.py_
~~~BASH
usage: opal.py [-h] -g GOLD_STANDARD_FILE [-l LABELS] -o OUTPUT_DIR [-r]
               profiles_files [profiles_files ...]

Compute all metrics for one or more taxonomic profiles

positional arguments:
  profiles_files        Files of profiles

optional arguments:
  -h, --help            show this help message and exit
  -g GOLD_STANDARD_FILE, --gold_standard_file GOLD_STANDARD_FILE
                        Gold standard file
  -l LABELS, --labels LABELS
                        Comma-separated profiles names
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Directory to write the results to
  -r, --by_rank         Create a results file per rank
~~~
**Example:**
~~~BASH
python3 opal.py -g data/goldstandard_low_1.bin \
data/cranky_wozniak_13 \
data/grave_wright_13 \
data/furious_elion_13 \
data/focused_archimedes_13 \
data/evil_darwin_13 \
data/agitated_blackwell_7 \
data/jolly_pasteur_3 \
-l "TIPP, Quikr, MP2.0, MetaPhyler, mOTU, CLARK, FOCUS" \
-o output_dir
~~~
**Output:**
Directory _output_dir_ will contain:
* a .tsv file for each profile (CLARK.tsv, FOCUS.tsv, MetaPhyler.tsv, mOTU.tsv, MP2.0.tsv, Quikr.tsv, and TIPP.tsv)
* spider_plot.pdf
* spider_plot_recall_precision.pdf
* plot_shannon.pdf

__Note 1__: spider plots will only be generated if at least 3 profiles are provided, so that the plots can form a triangle.

__Note 2__: _output_dir_ will will also contain a .tsv file per taxonomic rank if option -r is used.

# Developer Guide

We are using [tox]((https://tox.readthedocs.io/en/latest/)) for project automation.

### Tests

If you want to run tests, just type the following in the project's root directory:

~~~BASH
tox
~~~
