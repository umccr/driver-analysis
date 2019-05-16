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

...work in progress

[OncodriveFML](http://bbglab.irbbarcelona.org/oncodrivefml/home), see paper by [Mularoni et al., 2016](https://www.ncbi.nlm.nih.gov/pubmed/27311963) for details

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
