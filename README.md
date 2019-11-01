# Driver analysis

Workflow for **selection analysis** and **cancer driver discovery** using the following methods:

* **[dndscv](https://github.com/im3sanger/dndscv)** (see paper by [Martincorena et al., 2017](https://www.ncbi.nlm.nih.gov/pubmed/29056346) for details)
* **[OncodriveClust](https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html#92_detecting_cancer_driver_genes_based_on_positional_clustering)** (see paper by [Tamborero *et al*, 2013](https://www.ncbi.nlm.nih.gov/pubmed/23884480) for details)
* **[OncodriveFML](http://bbglab.irbbarcelona.org/oncodrivefml/home)** (see paper by [Mularoni *et al*, 2013](https://www.ncbi.nlm.nih.gov/pubmed/27311963) for details)
* **[Cancer Genome Interpreter](https://www.cancergenomeinterpreter.org/rest_api)** (CGI) (see paper by [Tamborero et al., 2018](https://www.ncbi.nlm.nih.gov/pubmed/29592813) for details for details)
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

For more details about [installation](OncodriveFML.md#installation), [data download](OncodriveFML.md#data-download), [configuration](OncodriveFML.md#Before-running-the-analysis), [usage](OncodriveFML.md#running-the-analysis) and [output](OncodriveFML.md#output) see documentaion in [OncodriveFML.md](OncodriveFML.md).

<br>

### Cancer Genome Interpreter

[Cancer Genome Interpreter](https://www.cancergenomeinterpreter.org) (CGI), see paper by [Tamborero et al., 2018](https://www.ncbi.nlm.nih.gov/pubmed/29592813) for details.

This tool is freely available through an API or a web interface at [http://www.cancergenomeinterpreter.org](http://www.cancergenomeinterpreter.org). 

The CGI resource can also be accessed programmatically by an API created via REST. Only registered users can make use of the API, since a token is needed for any communication between the end user and the REST API. Further details can be found at [https://www.cancergenomeinterpreter.org/api/v1](https://www.cancergenomeinterpreter.org/rest_api) and information about acceptable input format can be found at [https://www.cancergenomeinterpreter.org/formats](https://www.cancergenomeinterpreter.org/formats).

For more details about the [usage](OncodriveFML.md#running-the-analysis) and [output](OncodriveFML.md#output) see documentaion in [CGI.md](CGI.md).

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
--samples_id_cols | The name(s) of MAF file(s) column containing samples' IDs. One column name is expected for a single file, and each separated by comma. The defualt samples' ID column is `Tumor_Sample_Barcode` | No
--dnds_p | dN/dS method p-value threshold for reporting significant genes (defualt is `0.05`) | No
--dnds_q | dN/dS method q-value threshold for reporting significant genes (defualt is `1`) | No
--oncodriveclust_fdr | OncodriveClust method false discovery rate (FDR) threshold for reporting significant genes (defualt is `0.5`) | No
--ratios_ci | Calculate per-gene confidence intervals for the dN/dS ratios (default is `FALSE`) | No
--hypermut_sample_cutoff | Mutations per gene to define ultra-hypermutator samples (these will be excluded; defualt is `250`) | No
--max_muts_per_gene | Maximum mutations per gene in same sample (remaining will be subsampled; defualt is `3`) | No
--genes_list | Location and name of a file listing genes of interest to be considered in the report (OPTIONAL). NOTE, this option is implemented only in [dNdS](#dnds) method and relevant for targeted sequencing studies | No
--genes_blacklist | Location and name of a file listing genes to be excluded (OPTIONAL). Header is not expected and the genes should be listed in separate lines | No
--samples_blacklist | Location and name of a file listing samples to be excluded (OPTIONAL). The ID of samples to be excluded are expected to be listed in column named `Tumor_Sample_Barcode`. Additional columns are allowed | No
--nonSyn_list | List of variant classifications to be considered as non-synonymous. Rest will be considered as silent variants. Default uses [Variant Classifications](http://asia.ensembl.org/Help/Glossary?id=535) with `High/Moderate variant consequences` | No
--oncodrivefml | Name of folder with results files from [OncodriveFML](https://oncodrivefml.readthedocs.io/en/latest/oncodriveFML.html) analysis (OPTIONAL) | No
--oncodrivefml_p | P-value threshold for reporting [OncodriveFML](https://oncodrivefml.readthedocs.io/en/latest/oncodriveFML.html) results (OPTIONAL; defualt is `0.1`) | No
--oncodrivefml_q | Q-value threshold for reporting [OncodriveFML](https://oncodrivefml.readthedocs.io/en/latest/oncodriveFML.html) results (OPTIONAL; defualt is `0.001`) | No
--oncodrivefml_conf | Directory and name of [OncodriveFML](https://oncodrivefml.readthedocs.io/en/latest/oncodriveFML.html) configuration file (OPTIONAL) | No
--cgi | Name of folder and the results files from [Cancer Genome Interpreter](https://www.cancergenomeinterpreter.org) (CGI) analysis | No
--remove_duplicated_variants | Remove repeated variants in a particuar sample, mapped to multiple transcripts of same gene? (defulat is `TRUE`). **Note**, option `TRUE` removes all repeated variants as duplicated entries. `FALSE` results in keeping all of them) | No
--out_folder | Output folder (defualt is `Driver_analysis_report`) | No
--hide_code_btn | Hide the *Code* button allowing to show/hide code chunks in the final HTML report. Available options are: `TRUE` (default) and `FALSE` | No
--ucsc_genome_assembly | Version of UCSC genome assembly to be used as a reference (defualt is `19`, other available option is `38`) | No

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
