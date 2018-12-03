# Examples: running Bioboxes of taxonomic profilers and assessing their results with OPAL

The following examples show how to run Bioboxes of taxonomic profilers on different datasets, tracking their runtimes and maximum main memory usages, and automatically assessing their results.

If you have already run a profiler and want to assess its results, you only need to run [opal.py](README.md#running-opalpy) and can probably skip these examples.

The assessed taxonomic profilers, in these examples, are:

| Taxonomic profiler | Biobox docker image |
|---|---|
| CommonKmers | stefanjanssen/dockerprofilingtools:commonkmers |
| FOCUS 0.31 adapted for CAMI | stefanjanssen/dockerprofilingtools:focus |
| Quikr | stefanjanssen/dockerprofilingtools:quickr |
| mOTU 1.1 | stefanjanssen/dockerprofilingtools:motu |
| Metaphlan 2.2.0 | stefanjanssen/dockerprofilingtools:metaphlan2 |
| Metaphyler 1.25 | stefanjanssen/dockerprofilingtools:metaphyler |
| TIPP 2.0.0 | stefanjanssen/dockerprofilingtools:tipp |

## Comparing taxonomic profilers on the CAMI I high complexity dataset

- Download the 5 samples of the CAMI I high complexity dataset from <https://data.cami-challenge.org/participate> and save all files in the same directory. Your directory should contain:

~~~BASH
RH_S001__insert_270.fq.gz
RH_S002__insert_270.fq.gz
RH_S003__insert_270.fq.gz
RH_S004__insert_270.fq.gz
RH_S005__insert_270.fq.gz
~~~

- Pull the Bioboxes of profilers:

~~~BASH
docker pull stefanjanssen/docker_profiling_tools:commonkmers
docker pull stefanjanssen/docker_profiling_tools:focus
docker pull stefanjanssen/docker_profiling_tools:metaphlan2
docker pull stefanjanssen/docker_profiling_tools:metaphyler
docker pull stefanjanssen/docker_profiling_tools:quickr
docker pull stefanjanssen/docker_profiling_tools:tipp
docker pull stefanjanssen/docker_profiling_tools:motu
~~~

- CommonKmers uses a database that is not stored inside its Biobox. Download it from <https://zenodo.org/record/1749272/files/CommonKmersData.tar.gz?download=1> (DOI: <http://doi.org/10.5281/zenodo.1749272>) and extract the files. Make sure to set the path to the files with option `--volume` for `opal_workflow.py`, as shown below.

~~~BASH
wget --content-disposition https://zenodo.org/record/1749272/files/CommonKmersData.tar.gz?download=1
tar -xzf CommonKmersData.tar.gz
~~~

- OPAL's tool to run Bioboxes of profilers, measure their run time and maximum memory usage, and automatically assess their results is `opal_workflow.py`. To run it, you also need the gold standard file [gs_cami_i_hc.profile](data/gs_cami_i_hc.profile) and the Biobox YAML file [biobox_cami_i_hc.yaml](data/biobox_cami_i_hc.yaml).

- Run `opal_workflow.py` as follows, modifying the options to match your system's paths.

~~~BASH
python3 ./opal_workflow.py \
stefanjanssen/docker_profiling_tools:commonkmers \
stefanjanssen/docker_profiling_tools:focus \
stefanjanssen/docker_profiling_tools:metaphlan2 \
stefanjanssen/docker_profiling_tools:metaphyler \
stefanjanssen/docker_profiling_tools:quickr \
stefanjanssen/docker_profiling_tools:tipp \
stefanjanssen/docker_profiling_tools:motu \
--labels "CommonKmers, FOCUS, Metaphlan, MetaPhyler, Quikr, TIPP, mOTU" \
--input_dir /path/to/gzipped/fastq/files \
--output_dir /path/to/output_dir \
--yaml /path/to/biobox_cami_i_hc.yaml \
--volume /path/to/CommonKmersData:/exchange/db:ro \
--gold_standard_file data/gs_cami_i_hc.profile \
--plot_abundances \
--desc "1st CAMI Challenge Dataset 3 CAMI high"
~~~

The output directory, `output_dir` in this example, will be created if does not exist. It will contain the predictions of all profilers and OPAL's assessments.

## Comparing taxonomic profilers on the CAMI II mouse gut dataset

- Download the 64 short-read samples of the CAMI II mouse gut dataset from <https://data.cami-challenge.org/participate>. The files have the same name, but should be located in different sub-directories of the same root directory:

~~~BASH
2017.12.29_11.37.26_sample_0/reads/anonymous_reads.fq.gz
2017.12.29_11.37.26_sample_1/reads/anonymous_reads.fq.gz
2017.12.29_11.37.26_sample_2/reads/anonymous_reads.fq.gz
...
2017.12.29_11.37.26_sample_63/reads/anonymous_reads.fq.gz
~~~

- To run `opal_workflow.py`, you also need the gold standard file [gs_cami_i_hc.profile](data/gs_cami_i_hc.profile) and the Biobox YAML file [biobox_cami_ii_mg.yaml](data/biobox_cami_ii_mg.yaml).

- Follow and adapt the other steps given above.

## Comparing taxonomic profilers on the Human Microbiome Project Mock Community dataset

- Download the FASTQ file of the staggered sample (accession SRX055381) from NCBI SRA (<https://www.ncbi.nlm.nih.gov/sra>) and compress it using gzip. You should have file:

~~~BASH
SRR172903.fastq.gz
~~~

- To run `opal_workflow.py`, you also need the gold standard file [gs_hmp_mc.profile](data/gs_hmp_mc.profile) and the Biobox YAML file [biobox_hmp_mc.yaml](data/biobox_hmp_mc.yaml).

- Follow and adapt the other steps given above.