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

The analysis are executed using `oncodrivefml` command followed by [paramters](#parameters) of interest, e.g.

```
conda activate driver-analysis

data="examples"

cd $data

# Run OncodriveFML using MAF from ICGC PACA-AU samples
oncodrivefml --input simple_somatic_mutation.open.PACA-AU.maf --elements /g/data3/gx8/extras/jmarzec/apps//path/to/oncodrivefml/example/data/example/cds.tsv.gz --sequencing wgs --output  ICGC_PACA-AU_oncodrivefml_analysis

# Run OncodriveFML using MAF from ICGC PACA-CA samples
oncodrivefml --input simple_somatic_mutation.open.PACA-CA.maf --elements /g/data3/gx8/extras/jmarzec/apps//path/to/oncodrivefml/example/data/example/cds.tsv.gz --sequencing wgs --output  ICGC_PACA-CA_oncodrivefml_analysis
```

<br />

### Output

[OncodriveFML](http://bbglab.irbbarcelona.org/oncodrivefml/home) generates 3 output files with the same name but different extension. The name given to the files is the same as the name of the mutations file (*simple_somatic_mutation.open.PACA-AU.maf* and *simple_somatic_mutation.open.PACA-CA.maf* in the [example](#running-the-analysis) above) followed by `-oncodrivefml` and the extension:

* `.tsv` - tabulated file with the analysis results
* `.png` - an image with the most significant genes labeled
* `.html` - HTML file with an interactive plot which can be used to search for specific genes


<br />
