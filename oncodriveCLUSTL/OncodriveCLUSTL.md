# Driver analysis - OncodriveCLUSTL

## Table of contents

<!-- vim-markdown-toc GFM -->
* [Description](#description)
* [Installation](#installation)
* [Data download](OncodriveFML.md#data-download)
* [Before running the analysis](OncodriveFML.md#Before-running-the-analysis)
* [Running the analysis](#running-the-analysis)
* [Output](#output)

<!-- vim-markdown-toc -->

<br>

### Description

[OncodriveCLUSTL](http://bbglab.irbbarcelona.org/oncodriveclustl/home) (see paper by [Mularoni et al., 2016](https://academic.oup.com/bioinformatics/article/35/22/4788/5522012?login=true) for details) is a new nucleotide sequence-based clustering algorithm to detect cancer drivers in genomic regions. OncodriveCLUSTL is based on a local background model derived from the nucleotide context mutational probabilities of the cohort under study. Our method is able to identify well-known cancer drivers in coding regions. It can be applied to **non-coding regions** and non-human data.

### Installation

To install follow the [steps](https://bitbucket.org/bbglab/oncodriveclustl/src/master/). On Gadi the following steps have worked:

```
$ conda create -n oncodriveclustl-3.7 python=3.7
$ pip install oncodrivefml
```

The following command will show you the help:

```
$ oncodriveclustl --help
```

<br>

### Data download

The [scores and reference data](#reference-data) will be automatically downloaded the first time  *[OncodriveFML](https://bitbucket.org/bbglab/oncodrivefml/src/master/)* is run, but to speed up the process it is better to first download it using `BgData` package management tool that is installed together with *OncodriveFML*.

#### Example installation and run script

```
/g/data3/gx8/extras/sehrishk/miniconda/envs/oncodriveclustl-3.7/bin/oncodriveclustl

/g/data/gx8/projects/Kanwal_pdac_atlas/driver_analysis/run_oncodriveclustl.sh
```

### Before running the analysis

#### Parameters

OncodriveCLUSTL only requires two main inputs, the mutations file and the annotations file. Foloowing are the required inputs:

Argument | Description
------------ | ------------
--input-file | Mutations file with 5 required columns *[CHROMOSOME, POSITION, REF, ALT and SAMPLE](https://bitbucket.org/bbglab/oncodriveclustl/src/master/)*. Of note, The tsv file can be prepared from the MAF file by extracting and renaming `"sampleID" "chr"      "pos"      "ref"      "mut"`  coloumns.
--regions-file | File listing genomic elements to analyse. It requires 4 columns *[CHROMOSOME, START, END, ELEMENT]*
--genome | Genome to use. Default is hg19
--output | Name of the output folder

<br /> 

The `regions-file` for hg38 coding sequence can be preapared using [generate_cds_hg38.py](./generate_cds_hg38.py)


### Running the analysis

The analysis are executed using `oncodriveclustl` command followed by [paramters](#parameters) of interest, e.g.

```
#!/bin/bash
#PBS -P gx8
#PBS -q normalbw
#PBS -l walltime=48:00:00
#PBS -l mem=240GB
#PBS -l ncpus=7
#PBS -l wd
#PBS -l storage=gdata/gx8

# this generates an error if set. This is used by GISTIC tool
unset LD_LIBRARY_PATH

/g/data3/gx8/extras/sehrishk/miniconda/envs/oncodriveclustl-3.7/bin/oncodriveclustl --input-file /g/data/gx8/projects/Kanwal_pdac_atlas/driver_analysis/pdac_oncodriveclustl_analysis/pdac_samples.tsv --regions-file /g/data/gx8/projects/Kanwal_pdac_atlas/driver_analysis/output_cds_element.tsv.gz --genome hg38 --output-directory /g/data/gx8/projects/Kanwal_pdac_atlas/driver_analysis/pdac_cds_update_oncodrivefmoncodriveclustl
```

Running with default simulating, smoothing and clustering OncodriveCLUSTL parameters is not recommended as these may not be optimal for our data.
Supplementary Methods has details on how to perform model selection for specific data.

Same for Signatures which by default are calculated as mutation frequencies: # mutated ref>alt k-mer counts / # total substitutions
Supplementary Methods has information on how to perform a more accurate signatures calculation.

<br />

### Output

OncodriveCLUSTL generates 3 output files.

* Elements results file (`elements_results.txt`). TSV file containing results of the analyzed elements
* Clusters results file (`clusters_results.tsv`). TSV file containing results of the clusters observed in the analyzed elements
* Log file (`results.log`). TXT file containing OncodriveCLUSTL's run information


<br />
