[![CircleCI](https://circleci.com/gh/CAMI-challenge/OPAL.svg?style=svg)](https://circleci.com/gh/CAMI-challenge/OPAL)

# OPAL - Profiling Assessment

Example page produced by OPAL: *https://cami-challenge.github.io/OPAL/*


# Requirements

* python &ge; 3.5
* numpy &ge; 1.13.0
* matplotlib &ge; 2.0.2
* dendropy &ge; 4.3.0
* pandas &ge; 0.20.3
* biom-format &ge; 2.1.6
* scikit-bio &ge; 0.5.1
* bokeh &ge; 0.12.9
* seaborn &ge; 0.8.1

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
pip3 install biom-format==2.1.6 --user
pip3 install scikit-bio==0.5.1 --user
~~~

## Input
OPAL uses at least two files:
1. A gold standard taxonomic profile
2. One or more taxonomic profiles to be assessed

Files must be in the [CAMI profiling Bioboxes format](https://github.com/bioboxes/rfc/tree/master/data-format) or in the [BIOM (Biological Observation Matrix) format](http://biom-format.org/). Program [_tsv2biom.py_](#running-tsv2biompy) allows to convert profiles from the former format to the latter.

**The BIOM format**

The BIOM format used by OPAL is a sparse matrix stored in a JSON or HDF5 file, with a column per sample and a row per taxonomy ID, storing the corresponding abundances. RANK, TAXPATH, and TAXPATHSN are stored as metadata of each row and have the same meaning as in the CAMI profiling Bioboxes format:
* RANK: taxonomic rank
* TAXPATH and TAXPATHSN: path from the root of the taxonomy to the respective current taxon, including the current taxon, separated by a `|`. TAXPATH and TAXPATHSN contain identifiers and plain names, respectively, of the taxonomies. For more details and examples, see [CAMI profiling Bioboxes format](https://github.com/bioboxes/rfc/tree/master/data-format).

## Computed metrics

* Unifrac error
* L1 norm error
* True positives, false positives, false negatives
* Precision
* Recall
* F1 score
* Jaccard index
* Shannon diversity and equitability indices
* Bray–Curtis distance

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
* results.html
* results.tsv
* subdirectory per_rank with a .tsv file per taxonomic rank
* subdirectory per_tool with a .tsv file per tool (CLARK.tsv, FOCUS.tsv, MetaPhyler.tsv, mOTU.tsv, MP2.0.tsv, Quikr.tsv, and TIPP.tsv)
* spider_plot.pdf
* spider_plot_recall_precision.pdf
* plot_shannon.pdf

__Note__: spider plots will only be generated if at least 3 profiles are provided, so that the plots can form a triangle.

## Running _tsv2biom.py_
~~~BASH
usage: tsv2biom.py [-h] -o OUTPUT_FILE [-j] files [files ...]

Convert profile in the CAMI Bioboxes format to BIOM

positional arguments:
  files                 Input file(s), one file per sample

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Output file
  -j, --json            Output in json (default: hdf5)
~~~
**Example:**
~~~BASH
python3 tsv2biom.py data/cranky_wozniak_13 -o output_dir/cranky_wozniak_13.biom
~~~

# Developer Guide

We are using [tox]((https://tox.readthedocs.io/en/latest/)) for project automation.

### Tests

If you want to run tests, just type the following in the project's root directory:

~~~BASH
tox
~~~
