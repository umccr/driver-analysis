# Driver analysis - OncodriveFML

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

[OncodriveFML](http://bbglab.irbbarcelona.org/oncodrivefml/home) (see paper by [Mularoni et al., 2016](https://www.ncbi.nlm.nih.gov/pubmed/27311963) for details) is a method designed to analyse the pattern of somatic mutations across tumours in both **coding and non-coding genomic regions** to identify signals of positive selection, and therefore, their involvement in tumorigenesis. The identification of protein coding genes, promoters, untranslated regions, intronic splice regions, and lncRNAs-containing driver mutations in several malignancies using [OncodriveFML](https://oncodrivefml.readthedocs.io/en/latest/oncodriveFML.html) is described by [Mularoni et al., 2016](https://www.ncbi.nlm.nih.gov/pubmed/27311963).

### Installation

To install [OncodriveFML](http://bbglab.irbbarcelona.org/oncodrivefml/home) follow the [steps](https://bitbucket.org/bbglab/oncodrivefml/src/master/) below. Preferably, install [OncodriveFML](http://bbglab.irbbarcelona.org/oncodrivefml/home) within `driver-analysis` Conda environment.

```
pip install oncodrivefml
```

<br /> 

### Data download

The [scores and reference data](#reference-data) will be automatically downloaded the first time  *[OncodriveFML](https://bitbucket.org/bbglab/oncodrivefml/src/master/)* is run, but to speed up the process it is better to first download it using `BgData` package management tool that is installed together with *OncodriveFML*.

#### Example data

```
cd /g/data3/gx8/extras/jmarzec/apps

mkdir oncodrivefml

cd oncodrivefml

wget https://bitbucket.org/bbglab/oncodrivefml/downloads/oncodrivefml-examples_v2.2.tar.gz --no-check-certificate

tar -xvzf oncodrivefml-examples_v2.2.tar.gz
```

#### Reference data

Download reference data using `BgData` tool

###### Note

As default `bgdata` downloads the data to home directory (`~/.bgdata`). Change this directory in case the space on home directory is limited (see [bgdata docs](https://bgdata.readthedocs.io/en/latest/configuration.html) for more details).

`vim ~/.config/bbglab/bgdatav2.conf`

and change  `local_repository = "~/.bgdata"` to `"/g/data3/gx8/extras/jmarzec/apps/oncodrivefml/.bgdata"`

* precomputed Combined Annotation Dependent Depletion ([CADD](https://cadd.gs.washington.edu/info)) scores (~17Gb)

```
bgdata get genomicscores/caddpack/1.0
```

* genome reference (3Gb)

```
bgdata get datasets/genomereference/hg19
```

* gene stops (16Mb) (not required)
```
bgdata get datasets/genestops/hg19
```

###### Note

* If the reference data is located in directory other than home directory (e.i. `~/.bgdata`) then one needs to specify the location of the `scores file` in `~/.config/bbglab/oncodrivefml_v2.conf` (see [Configuration file](#configuration-file) section)

`vim ~/.config/bbglab/oncodrivefml_v2.conf`

and change `file = "%(bgdata://genomicscores/caddpack/1.0)"` to `"/g/data3/gx8/extras/jmarzec/apps/oncodrivefml/.bgdata/genomicscores/caddpack/1.0-20170217"`

* Additionally, it may be necessary to run oncodrivefml with the default `file = '%(bgdata://genomicscores/caddpack/1.0)'` for the **first time** since this triggers some changes in the background that enable the tool to search for reference files in directory other than home directory (e.i. `~/.bgdata`)

<br /> 

### Before running the analysis

#### Parameters

The key [parameters](https://oncodrivefml.readthedocs.io/en/latest/workflow.html#the-command-line-interface) to run [OncodriveFML](http://bbglab.irbbarcelona.org/oncodrivefml/home) are described below

Argument | Description
------------ | ------------
--input | Mutations file with 5 required columns *[CHROMOSOME, POSITION, REF, ALT and SAMPLE](https://oncodrivefml.readthedocs.io/en/latest/files.html#input-file-format)*. Of note, the program recognises the MAF-specific column names as well, i.e. `Chromosome`, `Start_Position`, `Reference_Allele`, `Tumor_Seq_Allele1` and `Tumor_Sample_Barcode`. Additional first row `"#version [version number]"` is fine as well
--elements | File listing genomic elements to analyse. It requires 4 columns *[CHROMOSOME, START, END, ELEMENT](https://oncodrivefml.readthedocs.io/en/latest/files.html#regions-file-format)*
--sequencing | Type of sequencing. Available options: `wgs`, `wes` and `targeted`. `targeted` option should be used if the platform is unknown
--output | Name of the output folder
--samples-blacklist | List of samples to be removed when loading the input file
--no-indels | Defines whether to discard indels in the analysis. Do not use this parameter if indels should be included
<br /> 

#### Configuration file

The remaining settings can be set up in [configuration file](https://oncodrivefml.readthedocs.io/en/latest/configuration.html).

By default, OncodriveFML is prepared to analyse mutations using **hg19** reference genome. For other genomes, update the configuration file accordingly (see [this](https://oncodrivefml.readthedocs.io/en/latest/configuration.html) guidline).

The following settings are defined as default:

Section | Parameter | Default value | Description
------------ | ------------ | ------------ | ------------
[genome]  | build | hg19 | Reference genome
[signature] | method | complement | Use a 96 matrix with the signatures complemented
[score] | file | ~/.bgdata/genomicscores/caddpack | Path to score file (downloaded using `bg-data get genomicscores/caddpack/1.0` command)
[statistic] | method | amean | Mathematical method to use to compare observed and simulated values (arithmetic mean)
[statistic] | discard_mnp | False | Do not use/use multiple bases substitutions (MNP) in the analysis (include them)
[[indels]] | include | TRUE | Include indels in the analysis
[settings] | cores | (commented out) | Use all the available cores

<br />

###### Note

If the precomputed Combined Annotation Dependent Depletion ([CADD](https://cadd.gs.washington.edu/info)) scores file (~17Gb) was downloaded using `BgData` tool then one can change the `[score]` parameter to the location with`bgdata` folder (`"/g/data3/gx8/extras/jmarzec/apps/oncodrivefml/.bgdata/genomicscores/caddpack/1.0-20170217"`on Raijin (see [Reference data](#reference-data) section)

<br />

### Running the analysis

Conda `driver-summary` (see [driver analysis installation](https://github.com/umccr/driver-analysis#installation) section) needs to be activated frist.

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

###### Note

* As default a file with coding sequence [regions](https://oncodrivefml.readthedocs.io/en/latest/files.html#regions-file-format) (CDS) downloaded from [OncodriveFML](http://bbglab.irbbarcelona.org/oncodrivefml/home) website is used.

* In case of `ImportError: pycurl: libcurl link-time ssl backend (openssl) is different from compile-time ssl backend (none/other)` error message install `pycurl` using `pip` (from `python-sdk` [GitHub issues](https://github.com/transloadit/python-sdk/issues/4#issuecomment-418120668))

```
pip install pycurl==7.43.0 --global-option=build_ext --global-option="-L/usr/local/opt/openssl/lib" --global-option="-I/usr/local/opt/openssl/include"
```

* In case of `ImportError: Something is wrong with the numpy installation. While importing we detected an older version of numpy in ['...']. One method of fixing this is to repeatedly uninstall numpy until none is found, then reinstall this version.` error message uninstall `numpy` and then re-install it using `pip` (from `python-sdk` [GitHub issues](https://github.com/transloadit/python-sdk/issues/4#issuecomment-418120668))

```
pip uninstall numpy
pip install numpy
```

<br />

### Output

[OncodriveFML](http://bbglab.irbbarcelona.org/oncodrivefml/home) generates 3 output files with the same name but different extension. The name given to the files is the same as the name of the mutations file (*simple_somatic_mutation.open.PACA-AU.maf* and *simple_somatic_mutation.open.PACA-CA.maf* in the [example](#running-the-analysis) above) followed by `-oncodrivefml` and the extension:

* `.tsv` - tabulated file with the analysis results
* `.png` - an image with the most significant genes labeled
* `.html` - HTML file with an interactive plot which can be used to search for specific genes


<br />
