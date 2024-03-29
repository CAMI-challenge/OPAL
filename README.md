[![CircleCI](https://circleci.com/gh/CAMI-challenge/OPAL.svg?style=shield)](https://circleci.com/gh/CAMI-challenge/OPAL) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/cami-opal/README.html)

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
* [CAMI II mouse gut toy dataset (using parameter --filter 1)](https://cami-challenge.github.io/OPAL/cami_ii_mg_filter1/)
* [CAMI II mouse gut toy dataset (using parameters --filter 1 --normalize)](https://cami-challenge.github.io/OPAL/cami_ii_mg_filter1_normalized)
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

OPAL 1.0.12 has been tested with Python 3.10 and 3.11.

See [requirements.txt](requirements.txt) for all dependencies.

### Steps

You can [install OPAL using Docker](#running-opalpy-using-docker), [Bioconda](http://bioconda.github.io/recipes/cami-opal/README.html), or as follows.

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
pip3 install cami-opal
~~~

Make sure to add OPAL to your PATH:

~~~BASH
echo 'PATH=$PATH:${HOME}/.local/bin' >> ~/.bashrc
source ~/.bashrc
~~~

## Inputs
_**`Note: Support for the BIOM format has been dropped (temporarily) in OPAL 1.0.4 due to incompatibility with Python 3.7.*.`**_

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
usage: opal.py -g GOLD_STANDARD_FILE -o OUTPUT_DIR [-n] [-f FILTER] [-p] [-l LABELS] [-t TIME] [-m MEMORY] [-d DESC] [-r RANKS] [--metrics_plot_rel METRICS_PLOT_REL]
               [--metrics_plot_abs METRICS_PLOT_ABS] [--silent] [-v] [-h] [-b BRANCH_LENGTH_FUNCTION] [--normalized_unifrac]
               profiles_files [profiles_files ...]

OPAL: Open-community Profiling Assessment tooL

required arguments:
  profiles_files        Files of profiles
  -g GOLD_STANDARD_FILE, --gold_standard_file GOLD_STANDARD_FILE
                        Gold standard file
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Directory to write the results to

optional arguments:
  -n, --normalize       Normalize samples
  -f FILTER, --filter FILTER
                        Filter out the predictions with the smallest relative abundances summing up to [FILTER]% within a rank
  -p, --plot_abundances
                        Plot abundances in the gold standard (can take some minutes)
  -l LABELS, --labels LABELS
                        Comma-separated profiles names
  -t TIME, --time TIME  Comma-separated runtimes in hours
  -m MEMORY, --memory MEMORY
                        Comma-separated memory usages in gigabytes
  -d DESC, --desc DESC  Description for HTML page
  -r RANKS, --ranks RANKS
                        Highest and lowest taxonomic ranks to consider in performance rankings, comma-separated. Valid ranks: superkingdom, phylum, class, order, family, genus, species,
                        strain (default:superkingdom,species)
  --metrics_plot_rel METRICS_PLOT_REL
                        Metrics for spider plot of relative performances, first character, comma-separated. Valid metrics: w:weighted Unifrac, l:L1 norm, c:completeness, p:purity, f:false
                        positives, t:true positives (default: w,l,c,p,f)
  --metrics_plot_abs METRICS_PLOT_ABS
                        Metrics for spider plot of absolute performances, first character, comma-separated. Valid metrics: c:completeness, p:purity, b:Bray-Curtis (default: c,p)
  --silent              Silent mode
  -v, --version         show program's version number and exit
  -h, --help            Show this help message and exit

UniFrac arguments:
  -b BRANCH_LENGTH_FUNCTION, --branch_length_function BRANCH_LENGTH_FUNCTION
                        UniFrac tree branch length function (default: "lambda x: 1/x", where x=tree depth)
  --normalized_unifrac  Compute normalized version of weighted UniFrac by dividing by the theoretical max unweighted UniFrac
~~~
**Example:** To run the example, please download the files given in the [_data_](https://github.com/CAMI-challenge/OPAL/tree/master/data) directory.
~~~BASH
./opal.py -g data/goldstandard_low_1.bin \
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
docker run -v $(pwd):/host opal \
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
* Meyer, F., Bremges, A., Belmann, P., Janssen, S., McHardy, A.C., and Koslicki, D. **Assessing taxonomic metagenome profilers with OPAL.** *Genome Biology*, 20, 51 (2019). [https://doi.org/10.1186/s13059-019-1646-y](https://doi.org/10.1186/s13059-019-1646-y)

Part of OPAL's functionality was described in the CAMI manuscript. Thus please also cite:
* Sczyrba, A., Hofmann, P., Belmann, P. et al. **Critical Assessment of Metagenome Interpretation—a benchmark of metagenomics software.** Nat Methods 14, 1063–1071 (2017). [https://doi.org/10.1038/nmeth.4458](https://doi.org/10.1038/nmeth.4458)

or

* Meyer, F., Fritz, A., Deng, ZL. et al. **Critical Assessment of Metagenome Interpretation: the second round of challenges.** Nat Methods 19, 429–440 (2022). [https://doi.org/10.1038/s41592-022-01431-4](https://doi.org/10.1038/s41592-022-01431-4)
