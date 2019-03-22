[![CircleCI](https://circleci.com/gh/CAMI-challenge/OPAL.svg?style=svg)](https://circleci.com/gh/CAMI-challenge/OPAL)

# OPAL: Open-community Profiling Assessment tooL

Taxonomic metagenome profilers predict the presence and relative abundance of microorganisms from shotgun sequence samples of DNA isolated directly from a microbial community. Over the past years, there has been an explosive growth of software and algorithms for this task, resulting in a need for more systematic comparisons of these methods based on relevant performance criteria. OPAL implements commonly used performance metrics, including those of the first challenge of the Initiative for the [Critical Assessment of Metagenome Interpretation (CAMI)](http://cami-challenge.org), together with convenient visualizations.

**Computed metrics**

* Unifrac error
* L1 norm error
* True positives, false positives, false negatives
* Precision
* Recall
* F1 score
* Jaccard index
* Shannon diversity and equitability indices
* Bray–Curtis distance

**Example pages produced by OPAL**
* [CAMI I high complexity challenge dataset](https://cami-challenge.github.io/OPAL/cami_i_hc/)
* [CAMI II mouse gut toy dataset](https://cami-challenge.github.io/OPAL/cami_ii_mg/)
* [Human Microbiome Project Mock Community dataset](https://cami-challenge.github.io/OPAL/hmp_mc/)

**See also**
* [Assessments of profiling submissions to the 1st CAMI Challenge](https://cami-challenge.github.io/OPAL/cami_i_challenge_submissions/)

# User Guide

* [Installation](#installation)
* [Inputs](#inputs)
* [Running opal.py](#running-opalpy)
* [Running opal.py using Docker](#running-opalpy-using-docker)
* [Running tsv2biom.py](#running-tsv2biompy)
* [Measuring runtime and maximum main memory usage](#measuring-runtime-and-maximum-main-memory-usage)
* [More examples](EXAMPLES.md)

## Installation

### Requirements

OPAL requires Python 3.5.

See [default.txt](requirements/default.txt) for all dependencies.

### Steps

You can run [OPAL using Docker (see below)](#running-opalpy-using-docker) or install it as follows.

Install pip if not already installed (tested on Linux Ubuntu 18.04):

~~~BASH
sudo apt install python3-pip
~~~
Should you receive the message `Unable to locate package python3-pip`, enter the following commands and repeat the previous step.

~~~BASH
sudo add-apt-repository universe
sudo apt update
~~~

Then run:

~~~BASH
pip3 install numpy==1.15.3
pip3 install cami-opal
~~~

Make sure to add OPAL to your PATH:

~~~BASH
echo 'PATH=$PATH:${HOME}/.local/bin' >> ~/.bashrc
source ~/.bashrc
~~~

## Inputs
OPAL uses at least two files:
1. A gold standard taxonomic profile
2. One or more taxonomic profiles to be assessed

Files must be in the [CAMI profiling Bioboxes format](https://github.com/bioboxes/rfc/tree/master/data-format) or in the [BIOM (Biological Observation Matrix) format](http://biom-format.org/). Program [_tsv2biom.py_](#running-tsv2biompy) allows to convert profiles from the former format to the latter.

**The BIOM format**

The BIOM format used by OPAL is a sparse matrix stored in a JSON or HDF5 file, with a column per sample and a row per taxonomy ID, storing the corresponding abundances. RANK, TAXPATH, and TAXPATHSN are stored as metadata of each row and have the same meaning as in the CAMI profiling Bioboxes format:
* RANK: taxonomic rank
* TAXPATH and TAXPATHSN: path from the root of the taxonomy to the respective current taxon, including the current taxon, separated by a `|`. TAXPATH and TAXPATHSN contain identifiers and plain names, respectively, of the taxonomies. For more details and examples, see [CAMI profiling Bioboxes format](https://github.com/bioboxes/rfc/tree/master/data-format).

## Running _opal.py_
~~~BASH
usage: opal.py -g GOLD_STANDARD_FILE -o OUTPUT_DIR [-n] [-p] [-l LABELS]
               [-t TIME] [-m MEMORY] [-d DESC] [--silent] [-v] [-h]
               profiles_files [profiles_files ...]

OPAL: Open-community Profiling Assessment tooL
                                                                                                                                            
required arguments:                                                                                                                            
  profiles_files        Files of profiles                                                                                                          
  -g GOLD_STANDARD_FILE, --gold_standard_file GOLD_STANDARD_FILE                                                                                     
                        Gold standard file                                                                                                              
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR                                                                                                                  
                        Directory to write the results to                                                                                                     
                                                                                                                                                                
optional arguments:                                                                                                                                               
  -n, --no_normalization                                                                                                                                            
                        Do not normalize samples                                                                                                                     
  -p, --plot_abundances                                                                                                                                                
                        Plot abundances in the gold standard (can take some
                        minutes)
  -l LABELS, --labels LABELS
                        Comma-separated profiles names
  -t TIME, --time TIME  Comma-separated runtimes in hours
  -m MEMORY, --memory MEMORY
                        Comma-separated memory usages in gigabytes
  -d DESC, --desc DESC  Description for HTML page
  --silent              Silent mode
  -v, --version         show program's version number and exit
  -h, --help            Show this help message and exit
~~~
**Example:** To run the example, please download the files given in the [_data_](https://github.com/CAMI-challenge/OPAL/tree/master/data) directory.
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

## Running _opal.py_ using Docker

Download or git-clone OPAL from GitHub. In OPAL's directory, build the Docker image with the command:

~~~BASH
docker build -t opal:latest .
~~~

_opal.py_ can then be run with the `docker run` command. Example:

~~~BASH
docker run -v /path/to/OPAL:/host opal:latest \
opal.py -g /host/data/goldstandard_low_1.bin \
/host/data/cranky_wozniak_13 \
/host/data/grave_wright_13 \
/host/data/furious_elion_13 \
/host/data/focused_archimedes_13 \
/host/data/evil_darwin_13 \
/host/data/agitated_blackwell_7 \
/host/data/jolly_pasteur_3 \
-l "TIPP, Quikr, MP2.0, MetaPhyler, mOTU, CLARK, FOCUS" \
-o /host/output_dir
~~~

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

## Measuring runtime and maximum main memory usage

To measure the runtime and maximum main memory usage of a taxonomic profiler using OPAL, it must be converted to a Biobox docker image. Several Bioboxes are already available on Docker Hub (see [Examples page](EXAMPLES.md)). 

To build your own Biobox, general instructions are available at <http://bioboxes.org/>. Most importantly, the Biobox of a profiler must satisfy specific input and output formats (see section [Inputs](#inputs) above). Helpful examples of scripts and Dockerfiles are available at <https://github.com/CAMI-challenge/docker_profiling_tools>.

OPAL's tools to measure runtime and maximum main memory usage are:

* `opal_stats.py:` Runs the Biobox of a taxonomic profiler and tracks its runtime and main memory usage.

* `opal_workflow.py:` Runs the Bioboxes of one of more taxonomic profilers, tracks their runtimes and main memory usages using `opal_stats.py`, and automatically assesses their results with `opal.py`.

See example usage of these tools in the [Examples page](EXAMPLES.md).

Runtimes and memory usages can also be manually provided to `opal.py` using options `--time` and `--memory`. They will then be incorporated in the results files and the HTML report.

## More examples

See [Examples page](EXAMPLES.md).

# Developer Guide

We are using [tox]((https://tox.readthedocs.io/en/latest/)) for project automation.

### Tests

If you want to run tests, just type the following in the project's root directory:

~~~BASH
tox
~~~

# Citation
Please cite:
* Fernando Meyer, Andreas Bremges, Peter Belmann, Stefan Janssen, Alice Carolyn McHardy, and David Koslicki (2019). **Assessing taxonomic metagenome profilers with OPAL.** *Genome Biology*, 20:51. doi:[10.1186/s13059-019-1646-y](https://doi.org/10.1186/s13059-019-1646-y)

Part of OPAL's functionality was described in the CAMI manuscript. Thus please also cite:
* Alexander Sczyrba, Peter Hofmann, Peter Belmann, et al. (2017). **Critical Assessment of Metagenome Interpretation—a benchmark of metagenomics software.** *Nature Methods*, 14, 11:1063–1071. doi:[10.1038/nmeth.4458](https://doi.org/10.1038/nmeth.4458)
