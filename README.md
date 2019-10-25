# Driver analysis

Workflow for **selection analysis** and **cancer driver discovery** using the following methods:

* **[dndscv](https://github.com/im3sanger/dndscv)** (see paper by [Martincorena et al., 2017](https://www.ncbi.nlm.nih.gov/pubmed/29056346) for details)
* **[OncodriveClust](https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html#92_detecting_cancer_driver_genes_based_on_positional_clustering)** (see paper by [Tamborero *et al*, 2013](https://www.ncbi.nlm.nih.gov/pubmed/23884480) for details)
* **[OncodriveFML](http://bbglab.irbbarcelona.org/oncodrivefml/home)** (see paper by [Mularoni *et al*, 2013](https://www.ncbi.nlm.nih.gov/pubmed/27311963) for details)
* **[Cancer Genome Interpreter](https://www.cancergenomeinterpreter.org/rest_api)** (CGI) (see paper by [Tamborero et al., 2018](https://www.ncbi.nlm.nih.gov/pubmed/29592813) for details for details, ...*work in progress*)
* **[MutSig](http://software.broadinstitute.org/cancer/cga/mutsig)** (see paper by [Lawrence et al., 2013](https://www.ncbi.nlm.nih.gov/pubmed/23770567) for details, ...*work in progress*)
* **[Hierarchical HotNet](https://github.com/raphael-group/hierarchical-hotnet)** (see preprint by [Reyna et al., 2017](https://www.ncbi.nlm.nih.gov/pubmed/30423088) for details, it is improved version of [HotNet2](https://github.com/raphael-group/hotnet2), ...*work in progress*)
* **[CHASMplus](https://github.com/KarchinLab/CHASMplus)** (see paper by [Tokheim and Karchin., preprint](https://www.biorxiv.org/content/10.1101/313296v4) for details, ...*work in progress*)

The results from individual tools are **summarised** and **visualised** using *[driverAnalysis.R](./scripts/driverAnalysis.R)* script. 


## Table of contents

<!-- vim-markdown-toc GFM -->
* [Installation](#installation)
* [Driver analysis tools](#driver-analysis-tools)
  * [dNdS](#dnds)
  * [OncodriveClust](#oncodriveclust)
  * [OncodriveFML](#oncodrivefml)
  * [Cancer Genome Interpreter](#cancer-genome-interpreter)
  * [MutSig](#mutsig)
  * [Hierarchical HotNet](#hierarchical-hotnet)
  * [CHASMplus](#chasmplus)
* [Driver analysis summary](#driver-analysis-summary)
    * [Usage](#usage)
    * [Arguments](#arguments)
    * [Examples](#examples)

<!-- vim-markdown-toc -->

<br>

### Installation

To summarise results from the individual [driver analysis tools](#driver-analysis-tools) run [Driver analysis summary](#driver-analysis-summary) R script. In order to create *conda* environment and install required packages run the [environment.yaml](envm/environment.yaml) file. The `-p` flag should point to the *miniconda* installation path. For instance, to create `driver-analysis` environment using *miniconda* installed in `/miniconda` directory run the following command:

```
conda env create -p /miniconda/envs/driver-analysis --file envm/environment.yaml
```

Activate created `driver-analysis` *conda* environment before running the pipeline

```
conda activate driver-analysis
```

###### Note

[OncodriveFML](#oncodrivefml) has to be installed separately. The instruction for installing this tool is provided [hehe](#oncodrivefml).

## Driver analysis tools

### dNdS

The ***[dnds](https://github.com/im3sanger/dndscv)*** (including *dndscv* and *dNdSloc* models, see paper by [Martincorena et al., 2017](https://www.ncbi.nlm.nih.gov/pubmed/29056346) for details) method is implemented in *[driverAnalysis.R](./scripts/driverAnalysis.R)* script described in [Driver analysis summary](#driver-analysis-summary) section.

<br>

### OncodriveClust

The **[OncodriveClust](https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html#92_detecting_cancer_driver_genes_based_on_positional_clustering)** (see paper by [Tamborero *et al*, 2013](https://www.ncbi.nlm.nih.gov/pubmed/23884480) for details) method is implemented in *[driverAnalysis.R](./scripts/driverAnalysis.R)* script described in [Driver analysis summary](#driver-analysis-summary) section.

<br>

### OncodriveFML

[OncodriveFML](http://bbglab.irbbarcelona.org/oncodrivefml/home) (see paper by [Mularoni et al., 2016](https://www.ncbi.nlm.nih.gov/pubmed/27311963) for details) is a method designed to analyse the pattern of somatic mutations across tumours in both **coding and non-coding genomic regions** to identify signals of positive selection, and therefore, their involvement in tumorigenesis. The identification of protein coding genes, promoters, untranslated regions, intronic splice regions, and lncRNAs-containing driver mutations in several malignancies using [OncodriveFML](https://oncodrivefml.readthedocs.io/en/latest/oncodriveFML.html) is described by [Mularoni et al., 2016](https://www.ncbi.nlm.nih.gov/pubmed/27311963).

#### Installation

To install [OncodriveFML](http://bbglab.irbbarcelona.org/oncodrivefml/home) follow the [steps](https://bitbucket.org/bbglab/oncodrivefml/src/master/) below. Preferably, install [OncodriveFML](http://bbglab.irbbarcelona.org/oncodrivefml/home) within `driver-analysis` Conda environment.

```
pip install oncodrivefml
```

<br /> 

#### Data download

The [scores and reference data](#reference-data) will be automatically downloaded the first time  *[OncodriveFML](https://bitbucket.org/bbglab/oncodrivefml/src/master/)* is run, but to speed up the process it is better to first download it using `BgData` package management tool that is installed together with *OncodriveFML*.

##### Example data

```
cd /g/data3/gx8/extras/jmarzec/apps

mkdir oncodrivefml

cd oncodrivefml

wget https://bitbucket.org/bbglab/oncodrivefml/downloads/oncodrivefml-examples_v2.2.tar.gz --no-check-certificate

tar -xvzf oncodrivefml-examples_v2.2.tar.gz
```

##### Reference data

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

#### Before running the analysis

##### Parameters

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

##### Configuration file

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

#### Running the analysis

Conda `driver-summary` (see [driver analysis installation](https://github.com/umccr/driver-analysis#installation) section) needs to be activated frist.

The analysis are executed using `oncodrivefml` command followed by [paramters](#parameters) of interest, e.g.

```
conda activate driver-analysis

data="examples"

cd $data

# Run OncodriveFML using MAF from ICGC PACA-AU samples
oncodrivefml --input simple_somatic_mutation.open.PACA-AU.maf --elements /path/to/oncodrivefml/example/data/cds.tsv.gz --sequencing wgs --output  ICGC_PACA-AU_oncodrivefml_analysis

# Run OncodriveFML using MAF from ICGC PACA-CA samples
oncodrivefml --input simple_somatic_mutation.open.PACA-CA.maf --elements /g/data3/gx8/extras/jmarzec/apps//path/to/oncodrivefml/example/data/example/cds.tsv.gz --sequencing wgs --output  ICGC_PACA-CA_oncodrivefml_analysis
```

###### Note

* As default a file with coding sequence [regions](https://oncodrivefml.readthedocs.io/en/latest/files.html#regions-file-format) (CDS) downloaded from [OncodriveFML](http://bbglab.irbbarcelona.org/oncodrivefml/home) website is used.

* In case of `ImportError: pycurl: libcurl link-time ssl backend (openssl) is different from compile-time ssl backend (none/other)` error message install `pycurl` using `pip` (from `python-sdk` [GitHub issues](https://github.com/transloadit/python-sdk/issues/4#issuecomment-418120668))

```
pip install pycurl==7.43.0 --global-option=build_ext --global-option="-L/usr/local/opt/openssl/lib" --global-option="-I/usr/local/opt/openssl/include"
```

<br />

#### Output

[OncodriveFML](http://bbglab.irbbarcelona.org/oncodrivefml/home) generates 3 output files with the same name but different extension. The name given to the files is the same as the name of the mutations file (*simple_somatic_mutation.open.PACA-AU.maf* and *simple_somatic_mutation.open.PACA-CA.maf* in the [example](#running-the-analysis) above) followed by `-oncodrivefml` and the extension:

* `.tsv` - tabulated file with the analysis results
* `.png` - an image with the most significant genes labeled
* `.html` - HTML file with an interactive plot which can be used to search for specific genes

<br>

### Cancer Genome Interpreter

...work in progress

[Cancer Genome Interpreter](http://software.broadinstitute.org/cancer/cga/mutsig) (CGI), see paper by [Tamborero et al., 2018](https://www.ncbi.nlm.nih.gov/pubmed/29592813) for details.

This tool is freely available through an API or a web interface at [http://www.cancergenomeinterpreter.org](http://www.cancergenomeinterpreter.org). 

The CGI resource can also be accessed programmatically by an API created via REST. Only registered users can make use of the API, since a token is needed for any communication between the end user and the REST API. Further details can be found at [https://www.cancergenomeinterpreter.org/api/v1](https://www.cancergenomeinterpreter.org/rest_api).

<br>

### MutSig

...work in progress

[MutSig](http://software.broadinstitute.org/cancer/cga/mutsig), see paper by [Lawrence et al., 2013](https://www.ncbi.nlm.nih.gov/pubmed/23770567) for details.

<br>

### Hierarchical HotNet

...work in progress

[Hierarchical HotNet](https://github.com/raphael-group/hierarchical-hotnet), see paper by [Reyna et al., 2015](https://www.ncbi.nlm.nih.gov/pubmed/30423088) for details.

<br>

### CHASMplus

...work in progress

[CHASMplus](https://github.com/KarchinLab/CHASMplus), see paper by [Tokheim and Karchin., preprint](https://www.biorxiv.org/content/10.1101/313296v4) for details.

<br>

## Driver analysis summary

To summarise and visualise results derived using various methods run the *[driverAnalysis.R](./scripts/driverAnalysis.R)* script. This script catches the arguments from the command line and passes them to the *[driverAnalysis.Rmd](./scripts/driverAnalysis.Rmd)* script to produce the html report, generate set of plots and tables summarising driver analysis results from user-defined tools.

This script also performs selection analysis and cancer driver discovery using ***[dnds](https://github.com/im3sanger/dndscv)*** (including *dndscv* and *dNdSloc* models, [Martincorena et al., 2017](https://www.ncbi.nlm.nih.gov/pubmed/29056346)) and ***[OncodriveClust](https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html#92_detecting_cancer_driver_genes_based_on_positional_clustering)*** ([Tamborero *et al*, 2013](https://www.ncbi.nlm.nih.gov/pubmed/23884480)) methods. The input mutation data is expected to be in [Mutation Annotation Format](https://software.broadinstitute.org/software/igv/MutationAnnotationFormat) (MAF) (see [MAF-summary](https://github.com/umccr/MAF-summary) repo for more detials), which is processed using *[maftools](https://www.bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html)* R package ([bioRxiv](http://dx.doi.org/10.1101/052662), [GitHub](https://github.com/PoisonAlien/maftools)).

###### Note

Only non-synonymous variants with high/moderate variant consequences, including *frame shift deletions*, *frame shift insertions*, *splice site mutations*, *translation start site mutations* ,*nonsense mutation*, *nonstop mutations*, *in-frame deletion*, *in-frame insertions* and *missense mutation*, are reported (silent variants are ignored). One can manually define variant classifications to be considered as non-synonymous using `--nonSyn_list` parameter.

While *dN/dS* and *OncodriveClust* methods use *[maftools](https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html)* object limited to non-synonymous variants ***OncodriveFML*** algorithm is able to analyse the pattern of somatic mutations across tumours in both **coding** and **non-coding** genomic regions to identify signals of positive selection. For that reason, certain information (e.g. in *Mutation maps* section) about variants not classified as non-synonymous will not be available in the summary report.


### Usage

To run the pipeline execure the *[driverAnalysis.R](./scripts/driverAnalysis.R)* script. This script catches the arguments from the command line and passes them to the *[driverAnalysis.Rmd](./scripts/driverAnalysis.Rmd)* script to produce the interactive HTML report.

#### Arguments

Argument | Description | Required
------------ | ------------ | ------------
--maf_dir | Directory with *MAF* file(s) | **Yes**
--maf_files | List of *MAF* file(s) to be processed. Each file name is expected to be separated by comma | **Yes**
--datasets | Desired names of each dataset. The names are expected to be in the same order as provided *MAF* files | **Yes**
--dnds_q | dN/dS method q-value threshold for reporting significant genes (defualt is `0.1`) | No
--oncodriveclust_fdr | OncodriveClust method false discovery rate (FDR) threshold for reporting significant genes (defualt is `0.5`) | No
--ratios_ci | Calculate per-gene confidence intervals for the dN/dS ratios (default is `FALSE`) | No
--hypermut_sample_cutoff | Mutations per gene to define ultra-hypermutator samples (these will be excluded; defualt is `250`) | No
--max_muts_per_gene | Maximum mutations per gene in same sample (remaining will be subsampled; defualt is `3`) | No
--ucsc_genome_assembly | Version of UCSC genome assembly to be used as a reference (defualt is `19`) | No
--genes_list | Location and name of a file listing genes of interest to be considered in the report (OPTIONAL). NOTE, this option is implemented only in [dNdS](#dnds) method and relevant for targeted sequencing studies | No
--genes_blacklist | Location and name of a file listing genes to be excluded (OPTIONAL). Header is not expected and the genes should be listed in separate lines | No
--samples_blacklist | Location and name of a file listing samples to be excluded (OPTIONAL). The ID of samples to be excluded are expected to be listed in column named `Tumor_Sample_Barcode`. Additional columns are allowed | No
--nonSyn_list | List of variant classifications to be considered as non-synonymous. Rest will be considered as silent variants. Default uses [Variant Classifications](http://asia.ensembl.org/Help/Glossary?id=535) with `High/Moderate variant consequences` | No
--oncodrivefml | Name of folder with results files from [OncodriveFML](https://oncodrivefml.readthedocs.io/en/latest/oncodriveFML.html) analysis (OPTIONAL) | No
--oncodrivefml_p | P-value threshold for reporting [OncodriveFML](https://oncodrivefml.readthedocs.io/en/latest/oncodriveFML.html) results (OPTIONAL; defualt is `0.1`) | No
--oncodrivefml_q | Q-value threshold for reporting [OncodriveFML](https://oncodrivefml.readthedocs.io/en/latest/oncodriveFML.html) results (OPTIONAL; defualt is `0.001`) | No
--oncodrivefml_conf | Directory and name of [OncodriveFML](https://oncodrivefml.readthedocs.io/en/latest/oncodriveFML.html) configuration file (OPTIONAL) | No
--remove_duplicated_variants | Remove repeated variants in a particuar sample, mapped to multiple transcripts of same gene? (defulat is `TRUE`). **Note**, option `TRUE` removes all repeated variants as duplicated entries. `FALSE` results in keeping all of them) | No
--out_folder | Output folder (defualt is `Driver_analysis_report`) | No
--hide_code_btn | Hide the *Code* button allowing to show/hide code chunks in the final HTML report. Available options are: `TRUE` (default) and `FALSE` | No

<br />

**Packages**: required packages are listed in [environment.yaml](envm/environment.yaml) file.

#### Examples 

Conda `driver-summary` (see [driver analysis installation](https://github.com/umccr/driver-analysis#installation) section) needs to be activated first.

Below are command line use examples for generating *Driver Analyses Summary* report using:

```
conda activate driver-analysis

oncodrivefml_conf=~/.config/bbglab/oncodrivefml_v2.conf

data="examples"

cd scripts

Rscript driverAnalysis.R --maf_dir $data --maf_files simple_somatic_mutation.open.PACA-AU.maf,simple_somatic_mutation.open.PACA-CA.maf --datasets ICGG_PACA-AU,ICGG_PACA-CA --dnds_q 0.1 --ratios_ci FALSE --hypermut_sample_cutoff 200 --max_muts_per_gene 3 --ucsc_genome_assembly 19 --oncodrivefml $data/ICGG_PACA-AU_oncodrivefml_analysis/simple_somatic_mutation.open.PACA-AU-oncodrivefml,$data/ICGG_PACA-CA_oncodrivefml_analysis/simple_somatic_mutation.open.PACA-AU-oncodrivefml --oncodrivefml_conf $oncodrivefml_conf --out_folder ICGC_PACA_analysis_report
```

<br />

This will generate *ICGC_PACA_analysis_report.html* report with summary tables and plots within `ICGC_PACA_analysis_report_report` folder, as well as `results` folder with intermediate files, including plots and tables that are presented in the report.
