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

The analysis are executed using `oncodriveclustl` command followed by [paramters](#parameters) of interest. Running with default simulating, smoothing and clustering OncodriveCLUSTL parameters is not recommended as these may not be optimal for our data.
Supplementary Methods has details on how to perform model selection for specific data.

Same for Signatures which by default are calculated as mutation frequencies: # mutated ref>alt k-mer counts / # total substitutions
Supplementary Methods has information on how to perform a more accurate signatures calculation.

Optimal parameters for OncoDriveCLUSTL were found empirically based on the methodolgy described in the original publication (Arnedo-Pac et al., 2019). For each combination of parameters, OncoDriveCLUSTL was run on PDAC dataset. All possible combinations of the following parameters were tested:

	- --smooth-window - [11 15 21 25 31 35 41 45]
	- --cluster-window - [11 15 21 25 31]
	- --simulation-window - [31 35 40]

The script used to generate combinations and run OncoDriveCLUSTL is below:

```
#!/bin/bash
#PBS -P gx8
#PBS -q normalbw
#PBS -l walltime=48:00:00
#PBS -l mem=240GB
#PBS -l ncpus=7
#PBS -l wd
#PBS -l storage=gdata/gx8

unset LD_LIBRARY_PATH
# Define the arrays of values for each window type
smoothing_windows=(11 15 21 25 31 35 41 45)
clustering_windows=(11 15 21 25 31)
simulation_windows=(31 35 40)

# Generate all combinations
window_combinations=()

for smoothing in "${smoothing_windows[@]}"; do
  for clustering in "${clustering_windows[@]}"; do
    for simulation in "${simulation_windows[@]}"; do
      window_combinations+=("$smoothing $clustering $simulation")
    done
  done
done

# Print or use combinations as needed
for windows in "${window_combinations[@]}"; do
  read -r smoothing clustering simulation <<< "$windows"
  echo "Smoothing Window: $smoothing, Clustering Window: $clustering, Simulation Window: $simulation"

  # Run oncodriveclustl command with the current window combination
  /g/data3/gx8/extras/sehrishk/miniconda/envs/oncodriveclustl-3.7/bin/oncodriveclustl --input-file /g/data/gx8/projects/Kanwal_pdac_atlas/driver_analysis/pdac_oncodriveclustl_analysis/pdac_cptac_samples_rm_randoms.tsv \
    --regions-file /g/data/gx8/projects/Kanwal_pdac_atlas/driver_analysis/output_cds_element.tsv.gz \
    --genome hg38 --simulation-window "$simulation" --smooth-window "$smoothing" --cluster-window "$clustering" --qqplot \
    --output-directory "/g/data/gx8/projects/Kanwal_pdac_atlas/driver_analysis/pdac_oncodriveclustl_results/sim${simulation}_smo${smoothing}_clus${clustering}"
done
```

### Output

OncodriveCLUSTL generates 3 output files.

* Elements results file (`elements_results.txt`). TSV file containing results of the analyzed elements
* Clusters results file (`clusters_results.tsv`). TSV file containing results of the clusters observed in the analyzed elements
* Log file (`results.log`). TXT file containing OncodriveCLUSTL's run information

<br />

**Note:**

In the output from these runs, CGC column was indicated as "Not Available". It's unclear from the documentation, how this can be avoided. After downloading the results locally, [check_cgc_genes.py](scripts/check_cgc_genes.py) can be used to update this column. The script sets CGC column to True, if GE in elements_results.txt exists in [COSMIC gene list](data/cosmic_grch38_v96.txt).

<br />

### Model Selection

OncodriveCLUSTL is an unsupervised method to identify clustering of mutations along the genomic sequence of GEs. However, the method resorts to three main hyperparameters that strongly determine its performace: the shape of identified clusters depends on the smoothing (i) and clustering (ii) windows; the simulation of mutations depends on a sampling or simulation window (iii) which defines the region where mutations are randomly distributed.

Model selection is done based on two criteria: 

	i) goodness of fit of observed p-values vs. the uniform distribution; 
	ii) enrichment of bona-fide known cancer elements in the ranking given by the method as an output.

**Goodness of fit of observed p-values vs. the uniform distribution**

Helper function `calculate_ks` for this specific step is in [here](./scripts/select_best_config.py). It takes the subset of observed p-values which aregreater than 0.1; second, it samples 1,000 of them to avoid sample size biases in our comparisons; third, calculates the KS statistic. It is recommended to select all configurations bearing an absolute value of the KS statistic up to 10% larger than the minimum KS statistic. This procedure lefts us with a set of most suitable configurations.

**Enrichment of bona-fide known cancer elements**

Helper function `evaluate_enrichment_method` for this specific step is in [select_best_config.py](./scripts/select_best_config.py). As indicated by the authors, for each configuration, this script calculates a CGC genes enrichment score for the top ranking genes (n=100 - default is 40). For each 1≤n≤100, it computes the proportion of CGC genes within the subset of genes with rank n or lower. Then, all these proportions are added up, albeit giving more weight to the terms arising from smaller sets. The configuration with highest enrichment score (E) was selected.

The output results are written to a txt file [calculate_ks.txt](./output/calculate_ks_enrich.txt). For PDAC data, the best configurations is:

```
Directory             sim40_smo41_clus25
KS Statistic                     0.30801
KS P-value                      0.000723
Number of P-values                    40
ENRICH_LOG                      4.789501
```

<br />

