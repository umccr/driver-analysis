################################################################################
#
#   File name: driverAnalysis.R
#
#   Authors: Jacek Marzec ( jacek.marzec@unimelb.edu.au )
#
#   University of Melbourne Centre for Cancer Research,
#   Victorian Comprehensive Cancer Centre
#   305 Grattan St, Melbourne, VIC 3000
#
################################################################################

################################################################################
#
#   Description: Script for selection analyses and cancer driver discovery results. This script catches the arguments from the command line and passes them to the driverAnalysis.Rmd script to produce the report, generate set of plots and tables.
#
#   Command line use example: Rscript driverAnalysis.R --maf_dir /data --maf_filessimple_somatic_mutation.open.PACA-AU.maf,PACA-CA.icgc.simple_somatic_mutation.maf --datasets ICGC-PACA-AU,ICGC-PACA-CA --dnds_q 0.1 --ratios_ci FALSE --hypermut_sample_cutoff 1000 --max_muts_per_gene 3 --ucsc_genome_assembly 19 --out_folder Driver_analysis_report
#
#   maf_dir:      Directory with MAF files
#   maf_files:    List of MAF files to be processed. Each file name is expected to be separated by comma
#   datasets:     Desired names of each dataset. The names are expected to be in the same order as provided MAF files and should be separated by comma
#   dnds_q:       dNdS method q-value threshold for reporting significant genes (defualt 0.1)
#   oncodriveclust_fdr:   OncodriveClust method false discovery rate (FDR) threshold for reporting significant genes (defualt 0.5)
#   ratios_ci:    Calculate per-gene confidence intervals for the dN/dS ratios (default FALSE)
#   hypermut_sample_cutoff:   Mutations per gene to define ultra-hypermutator samples (these will be excluded; defualt 1000)
#   max_muts_per_gene:   Maximum mutations per gene in same sample (remaining will be subsampled; defualt 3)
#   ucsc_genome_assembly:   Version of UCSC genome assembly to be used as a reference
#   out_folder:   Name for the output folder that will be created within the directory with MAF files. If no output folder is specified the results will be saved in folder "Driver_analysis_report"
#   genes_list (optional):  Location and name of a file listing genes of interest to be considered in the report. The genes are expected to be listed in first column
#   genes_list (optional):  Location and name of a file listing genes of interest to be considered in the report. The genes are expected to be listed in first column
#   genes_blacklist (optional):  Location and name of a file listing genes to be excluded. Header is not expected and the genes should be listed in separate lines
#   nonSyn_list (optional):   List of variant classifications to be considered as non-synonymous. Rest will be considered as silent variants
#   oncodrivefml (optional):   Name of folder and the results files from OncodriveFML analysis
#   oncodrivefml_p (optional):   P-value threshold for reporting OncodriveFML results. Defualt values is 0.1
#   oncodrivefml_q (optional):   Q-value threshold for reporting OncodriveFML results. Defualt values is 0.001
#   oncodrivefml_conf (optional):   Directory and name of OncodriveFML configuration file
#	  remove_duplicated_variants (optional):		Remove repeated variants in a particuar sample, mapped to multiple transcripts of same gene? Defulat value is "FALSE"
#   hide_code_btn : Hide the "Code" button allowing to show/hide code chunks in the final HTML report. Available options are: "TRUE" (default) and "FALSE"
#
################################################################################

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()


#===============================================================================
#    Load libraries
#===============================================================================

suppressMessages(library(optparse))


#===============================================================================
#    Catching the arguments
#===============================================================================
option_list <- list(
  make_option(c("-d", "--maf_dir"), action="store", default=NA, type='character',
              help="Directory with MAF files"),
  make_option(c("-m", "--maf_files"), action="store", default=NA, type='character',
              help="List of MAF files to be processed"),
  make_option(c("-c", "--datasets"), action="store", default=NA, type='character',
              help="Desired names of each dataset"),
  make_option(c("-q", "--dnds_q"), action="store", default=0.1, type='double',
              help="dNdS method q-value threshold for reporting significant genes"),
  make_option(c("-k", "--oncodriveclust_fdr"), action="store", default=0.5, type='double',
              help="OncodriveClust method false discovery rate (FDR) threshold for reporting significant genes"),
  make_option(c("-r", "--ratios_ci"), action="store", default=FALSE, type='logical',
              help="Calculate per-gene confidence intervals for the dN/dS ratios"),
  make_option(c("-u", "--hypermut_sample_cutoff"), action="store", default=250, type='integer',
              help="Mutations per gene to define ultra-hypermutator samples"),
  make_option(c("-s", "--max_muts_per_gene"), action="store", default=3, type='integer',
              help="Maximum mutations per gene in same sample"),
  make_option(c("-g", "--ucsc_genome_assembly"), action="store", default=19, type='integer',
              help="Version of UCSC genome assembly to be used as a reference"),
  make_option(c("-o", "--out_folder"), action="store", default="Driver_analysis_report", type='character',
              help="Output directory"),
  make_option(c("-l", "--genes_list"), action="store", default=NA, type='character',
              help="Location and name of a file listing genes of interest to be considered in the report"),
  make_option(c("-b", "--genes_blacklist"), action="store", default=NA, type='character',
              help="Location and name of a file listing genes to be excluded"),
  make_option(c("-i", "--samples_blacklist"), action="store", default=NA, type='character',
              help="Location and name of a file listing samples to be excluded"),
  make_option(c("-n", "--nonSyn_list"), action="store", default=NA, type='character',
              help="List of variant classifications to be considered as non-synonymous"),
  make_option(c("-f", "--oncodrivefml"), action="store", default=NA, type='character',
              help="Name of folder and the results files from OncodriveFML analysis"),
  make_option(c("-e", "--oncodrivefml_p"), action="store", default=0.01, type='double',
              help="P-value threshold for reporting OncodriveFML results. Defualt values is 0.1"),
  make_option(c("-j", "--oncodrivefml_q"), action="store", default=0.1, type='double',
              help="Q-value threshold for reporting OncodriveFML results. Defualt values is 0.01"),
  make_option(c("-a", "--oncodrivefml_conf"), action="store", default=NA, type='character',
              help="Directory and name of OncodriveFML configuration file"),
  make_option(c("-v", "--remove_duplicated_variants"), action="store", default=TRUE, type='logical',
              help="Remove repeated variants in a particuar sample, mapped to multiple transcripts of same gene?"),
  make_option(c("-hc", "--hide_code_btn"), action="store", default=TRUE, type='logical',
              help="Hide the \"Code\" button allowing to show/hide code chunks in the final HTML report")
  )

opt <- parse_args(OptionParser(option_list=option_list))

##### Collect MAF files and correspondiong datasets names
opt$maf_files <- gsub("\\s","", opt$maf_files)
opt$datasets <- gsub("\\s","", opt$datasets)

##### Read in argument from command line and check if all were provide by the user
if (is.na(opt$maf_dir) || is.na(opt$maf_files) || is.na(opt$datasets) ) {

  cat("\nPlease type in required arguments!\n\n")
  cat("\ncommand example:\n\nRscript driverAnalysis.R --maf_dir /data --maf_filessimple_somatic_mutation.open.PACA-AU.maf,PACA-CA.icgc.simple_somatic_mutation.maf --datasets ICGC-PACA-AU,ICGC-PACA-CA --dnds_q 0.1 --ratios_ci FALSE --hypermut_sample_cutoff 1000 --max_muts_per_gene 3 --ucsc_genome_assembly 19 --out_folder Driver_analysis_report\n\n")
  q()
  
} else if ( length(unlist(strsplit(opt$maf_files, split=',', fixed=TRUE))) != length(unlist(strsplit(opt$datasets, split=',', fixed=TRUE))) ) {

  cat("\nMake sure that the number of datasets names match the number of queried MAF files\n\n")
  q()
}

##### Pre-define list of variant classifications to be considered as non-synonymous. Rest will be considered as silent variants. Default uses Variant Classifications with High/Moderate variant consequences (http://asia.ensembl.org/Help/Glossary?id=535)
if ( is.na(opt$nonSyn_list) ) {
  opt$nonSyn_list<- c("Frame_Shift_Del","Frame_Shift_Ins","Splice_Site","Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation")
} else {
  opt$nonSyn_list <- unlist(strsplit(opt$nonSyn_list, split=',', fixed=TRUE))
}

##### Pass the user-defined argumentas to the driverAnalysis.R markdown script and run the analysis
rmarkdown::render(input = "driverAnalysis.Rmd", output_dir = paste(opt$maf_dir, opt$out_folder, "Report", sep = "/"), output_file = paste0(opt$out_folder, ".html"), params = list(maf_dir = opt$maf_dir, maf_files = opt$maf_files, datasets = opt$datasets, dnds_q = opt$dnds_q, oncodriveclust_fdr = opt$oncodriveclust_fdr, ratios_ci = opt$ratios_ci, hypermut_sample_cutoff = opt$hypermut_sample_cutoff, max_muts_per_gene = opt$max_muts_per_gene, ucsc_genome_assembly = opt$ucsc_genome_assembly, out_folder = opt$out_folder, genes_list = opt$genes_list, genes_blacklist = opt$genes_blacklist, samples_blacklist = opt$samples_blacklist, nonSyn_list = opt$nonSyn_list, oncodrivefml = opt$oncodrivefml, oncodrivefml_p = opt$oncodrivefml_p, oncodrivefml_q = opt$oncodrivefml_q, oncodrivefml_conf = opt$oncodrivefml_conf, remove_duplicated_variants = opt$remove_duplicated_variants, hide_code_btn = opt$hide_code_btn))
