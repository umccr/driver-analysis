# Driver analysis

Workflow for **selection analysis** and **cancer driver discovery** using the following methods:

* **[dndscv](https://github.com/im3sanger/dndscv)** (see paper by [Martincorena et al., 2017](https://www.ncbi.nlm.nih.gov/pubmed/29056346) for details)
* **[OncodriveClust](https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html#92_detecting_cancer_driver_genes_based_on_positional_clustering)** (see paper by [Tamborero *et al*, 2013](https://www.ncbi.nlm.nih.gov/pubmed/23884480) for details)
* **[OncodriveFML](http://bbglab.irbbarcelona.org/oncodrivefml/home)** (see paper by [Mularoni *et al*, 2013](https://www.ncbi.nlm.nih.gov/pubmed/27311963) for details, ...*work in progress*)
* **[MutSig](http://software.broadinstitute.org/cancer/cga/mutsig)** (see paper by [Lawrence et al., 2013](https://www.ncbi.nlm.nih.gov/pubmed/23770567) for details, ...*work in progress*)
* **[CHASMplus](https://github.com/KarchinLab/CHASMplus)** (see paper by [Tokheim and Karchin., preprint](https://www.biorxiv.org/content/10.1101/313296v4) for details, ...*work in progress*)
* **[Hierarchical HotNet](https://github.com/raphael-group/hierarchical-hotnet)** (see preprint by [Reyna et al., 2017](https://www.ncbi.nlm.nih.gov/pubmed/30423088) for details, it is improved version of [HotNet2](https://github.com/raphael-group/hotnet2), ...*work in progress*)

The results from individual tools are **summarised** and **visualised** using *[driverAnalysis.R](./scripts/driverAnalysis.R)* script. 


## Table of contents

<!-- vim-markdown-toc GFM -->
* [Driver analysis tools](#driver-analysis-tools)
  * [dNdS](#dnds)
  * [OncodriveClust](#oncodriveclust)
  * [OncodriveFML](#oncodrivefml)
  * [MutSig](#mutsig)
  * [CHASMplus](#chasmplus)
  * [Hierarchical HotNet](#hierarchical-hotnet)
* [Driver analysis summary](#driver-analysis-summary)

<!-- vim-markdown-toc -->
<br>

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

To install [OncodriveFML](http://bbglab.irbbarcelona.org/oncodrivefml/home) on Spartan follow the steps below (from [here](https://bitbucket.org/bbglab/oncodrivefml/src/master/))

```
module load python-dateutil

pip install oncodrivefml
```

<br /> 

#### Data download

##### Example data

```
cd $HOME/applications

mkdir oncodrivefml

cd oncodrivefml

wget https://bitbucket.org/bbglab/oncodrivefml/downloads/oncodrivefml-examples_v2.2.tar.gz --no-check-certificate

tar -xvzf oncodrivefml-examples_v2.2.tar.gz
```

##### Reference data

Download reference data using `BgData` tool

* precomputed CADD scores (~17Gb)

```
bg-data get genomicscores/caddpack/1.0
```

* genome reference (3Gb)

```
bg-data get datasets/genomereference/hg19
```

* gene stops (16Mb) (not required)
```
bg-data get datasets/genestops/hg19
```

<br /> 

#### Before running the analysis

##### Parameters

The key [parameters](https://oncodrivefml.readthedocs.io/en/latest/workflow.html#the-command-line-interface) to run [OncodriveFML](http://bbglab.irbbarcelona.org/oncodrivefml/home) are described below

Argument | Description
------------ | ------------
--input | Mutations file with 5 required columns *[CHROMOSOME, POSITION, REF, ALT and SAMPLE]*(https://oncodrivefml.readthedocs.io/en/latest/files.html#input-file-format). Of note, the program seems to recognise the MAF-specific column names as well, i.e. Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele1 and Tumor_Sample_Barcode. Additional first row "#version 2.4" is fine as well
--elements | File listing genomic elements to analyse. It requires 4 columns *[CHROMOSOME, START, END, ELEMENT](https://oncodrivefml.readthedocs.io/en/latest/files.html#regions-file-format)
--sequencing | Type of sequencing. Available options: "wgs", "wes" and "targeted". "targeted" option should be used if the platform is unknown
--output | Name of the output folder
--samples-blacklist | List of samples to be removed when loading the input file
-no-indels | Discard indels in the analysis. Do not use this parameter if indels should be included
<br /> 

##### Configuration file

The remaining setting can be set up in [configuration file](https://oncodrivefml.readthedocs.io/en/latest/configuration.html).

By default, OncodriveFML is prepared to analyse mutations using **hg19** reference genome. For other genomes, update the configuration file accordingly (see [this](https://oncodrivefml.readthedocs.io/en/latest/configuration.html) guidline).

The following settings are defined as default

Section | Parameter | Value | Description
------------ | ------------ | ------------ | ------------
genome  | build | hg19 | Reference genome
signature | method | complement | Use a 96 matrix with the signatures complemented
score | file | ~/.bgdata/genomicscores/caddpack | Path to score file (downloaded using `bg-data get genomicscores/caddpack/1.0` command)
statistic | method | amean | Mathematical method to use to compare observed and simulated values (arithmetic mean)
statistic | discard_mnp | False | Do not use/use multiple bases substitutions (MNP) in the analysis (include them)
indels | include | TRUE | Include indels in the analysis
settings | cores | (commented out) | Use all the available cores

<br />

#### Running the analysis

The analysis are executed using `oncodrivefml` command followed by paramters of interest, e.g.

```
cd /data

oncodrivefml --input simple_somatic_mutation.open.PACA-AU.maf --elements $HOME/applications/oncodrivefml/example/cds.tsv.gz --sequencing wgs --output  PACA-AU_oncodrivefml_analysis
```

Note, as default a file with of coding sequence [regions](https://oncodrivefml.readthedocs.io/en/latest/files.html#regions-file-format) (CDS) downloaded from [OncodriveFML](http://bbglab.irbbarcelona.org/oncodrivefml/home) website is used.

<br>

### MutSig

...work in progress

[MutSig](http://software.broadinstitute.org/cancer/cga/mutsig), see paper by [Lawrence et al., 2013](https://www.ncbi.nlm.nih.gov/pubmed/23770567) for details.

<br>

### CHASMplus

...work in progress

[CHASMplus](https://github.com/KarchinLab/CHASMplus), see paper by [Tokheim and Karchin., preprint](https://www.biorxiv.org/content/10.1101/313296v4) for details.

<br>

### Hierarchical HotNet

...work in progress

[Hierarchical HotNet](https://github.com/raphael-group/hierarchical-hotnet), see paper by [Reyna et al., 2015](https://www.ncbi.nlm.nih.gov/pubmed/30423088) for details.





<br>

## Driver analysis summary

To summarise and visualise results derived using various methods run the *[driverAnalysis.R](./scripts/driverAnalysis.R)* script. This script catches the arguments from the command line and passes them to the *[driverAnalysis.Rmd](./scripts/driverAnalysis.Rmd)* script to produce the html report, generate set of plots and tables summarising driver analysis results from user-defined tools.

This script also performs selection analysis and cancer driver discovery using ***[dnds](https://github.com/im3sanger/dndscv)*** (including *dndscv* and *dNdSloc* models, [Martincorena et al., 2017](https://www.ncbi.nlm.nih.gov/pubmed/29056346)) and ***[OncodriveClust](https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html#92_detecting_cancer_driver_genes_based_on_positional_clustering)*** ([Tamborero *et al*, 2013](https://www.ncbi.nlm.nih.gov/pubmed/23884480)) methods. The input mutation data is expected to be in [Mutation Annotation Format](https://software.broadinstitute.org/software/igv/MutationAnnotationFormat) (MAF) (see [MAF-summary](https://github.com/umccr/MAF-summary) repo for more detials), which is processed using *[maftools](https://www.bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html)* R package ([bioRxiv](http://dx.doi.org/10.1101/052662), [GitHub](https://github.com/PoisonAlien/maftools)).

**Script**: *[driverAnalysis.R](./scripts/driverAnalysis.R)*

Argument | Description
------------ | ------------
--maf_dir | Directory with *MAF* file(s)
--maf_files | List of *MAF* file(s) to be processed. Each file name is expected to be separated by comma
--datasets | Desired names of each dataset. The names are expected to be in the same order as provided *MAF* files
--q_value | Q-value threshold for reporting significant genes (defualt 0.1)
--ratios_ci | Calculate per-gene confidence intervals for the dN/dS ratios (default FALSE)
--hypermut_sample_cutoff | Mutations per gene to define ultra-hypermutator samples (these will be excluded; defualt 250)
--max_muts_per_gene | Maximum mutations per gene in same sample (remaining will be subsampled; defualt 3)
--ucsc_genome_assembly | Version of UCSC genome assembly to be used as a reference (defualt 19)
--out_folder | Output folder (defualt "Driver_analysis_report")
--genes_list | Location and name of a file listing genes of interest to be considered in the report (OPTIONAL)
--genes_blacklist | Location and name of a file listing genes to be excluded (OPTIONAL). Header is not expected and the genes should be listed in separate lines
--samples_blacklist | Location and name of a file listing samples to be excluded (OPTIONAL). The ID of samples to be excluded are expected to be listed in column named "Tumor_Sample_Barcode". Additional columns are allowed
--oncodrivefml | Name of folder with results files from [OncodriveFML](https://oncodrivefml.readthedocs.io/en/latest/oncodriveFML.html) analysis (OPTIONAL)
--oncodrivefml_p | P-value threshold for reporting [OncodriveFML](https://oncodrivefml.readthedocs.io/en/latest/oncodriveFML.html) results (OPTIONAL). Defualt values is 0.1
--oncodrivefml_q | Q-value threshold for reporting [OncodriveFML](https://oncodrivefml.readthedocs.io/en/latest/oncodriveFML.html) results (OPTIONAL). Defualt values is 0.001
--oncodrivefml_conf | Directory and name of [OncodriveFML](https://oncodrivefml.readthedocs.io/en/latest/oncodriveFML.html) configuration file (OPTIONAL)
<br />

**Packages**: *[maftools](https://www.bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html)*, *[optparse](https://cran.r-project.org/web/packages/optparse/optparse.pdf)*, *[knitr](https://cran.r-project.org/web/packages/knitr/knitr.pdf)*, *[DT](https://rstudio.github.io/DT/)*, *[ggplot2](https://cran.r-project.org/web/packages/ggplot2/ggplot2.pdf)*, *[dndscv](http://htmlpreview.github.io/?http://github.com/im3sanger/dndscv/blob/master/vignettes/dNdScv.html)*, *[UpSetR](https://cran.r-project.org/web/packages/UpSetR/README.html)*, *[stringr](https://cran.r-project.org/web/packages/stringr/vignettes/stringr.html)*, *[magick](https://cran.r-project.org/web/packages/magick/vignettes/intro.html)*

**Command line use example**:

```
Rscript driverAnalysis.R --maf_dir /data --maf_files simple_somatic_mutation.open.PACA-AU.maf,simple_somatic_mutation.open.PACA-CA.maf --datasets ICGC-PACA-AU,ICGC-PACA-CA --q_value 0.1 --ratios_ci FALSE --hypermut_sample_cutoff 3000 --max_muts_per_gene 3 --ucsc_genome_assembly 19 --out_folder Driver_analysis_report
```
<br>

This will generate *Driver_analysis_report.html* report with summary tables and plots within *Driver_analysis_report* folder.
