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
#   Command line use example: Rscript driverAnalysis.R --maf_dir /data --maf_filessimple_somatic_mutation.open.PACA-AU.maf,PACA-CA.icgc.simple_somatic_mutation.maf --datasets ICGC-PACA-AU,ICGC-PACA-CA --dnds_p 0.05 --ratios_ci FALSE --hypermut_sample_cutoff 3000 --max_muts_per_gene 3 --ucsc_genome_assembly 19 --out_folder Driver_analysis_report
#
#   maf_dir:      Directory with MAF files
#   maf_files:    List of MAF files to be processed. Each file name is expected to be separated by comma
#   datasets:     Desired names of each dataset. The names are expected to be in the same order as provided MAF files and should be separated by comma
#   samples_id_cols:  The name(s) of MAF file(s) column containing samples' IDs. One column name is expected for a single file, and each separated by comma. The defualt samples' ID column is "Tumor_Sample_Barcode"
#   dnds_p:       dNdS method p-value threshold for reporting significant genes (defualt 0.05)
#   dnds_q:       dNdS method q-value threshold for reporting significant genes (defualt 1)
#   activedriverwgs_p:       ActiveDriverWGS method p-value threshold for reporting significant genes (defualt 0.05)
#   activedriverwgs_fdr:       ActiveDriverWGS method FDR-value threshold for reporting significant genes (defualt 1)
#   activedriverwgs_cores:    Number of cores to be used for running ActiveDriverWGS method (defualt 1)
#   activedriverwgs_all_genes:     Run ActiveDriverWGS method for all genes (defualt FALSE where the analysed regions will be limited to genes listed in cancer genes (https://github.com/vladsaveliev/NGS_Utils/blob/master/ngs_utils/reference_data/key_genes/umccr_cancer_genes.2019-03-20.tsv))
#   oncodriveclust_fdr:   OncodriveClust method false discovery rate (FDR) threshold for reporting significant genes (defualt 0.5)
#   ratios_ci:    Calculate per-gene confidence intervals for the dN/dS ratios (default FALSE)
#   hypermut_sample_cutoff:   Mutations per gene to define ultra-hypermutator samples (these will be excluded; defualt 1000)
#   max_muts_per_gene:   Maximum mutations per gene in same sample (remaining will be subsampled; defualt 3)
#   genes_list (optional):  Location and name of a file listing genes of interest to be considered in the report. The genes are expected to be listed in first column
#   genes_blacklist (optional):  Location and name of a file listing genes to be excluded. Header is not expected and the genes should be listed in separate lines
#   nonSyn_list (optional):   List of variant classifications to be considered as non-synonymous. Rest will be considered as silent variants
#   oncodrivefml (optional):   Name of folder and the results files from OncodriveFML analysis
#   oncodrivefml_p (optional):   P-value threshold for reporting OncodriveFML results. Defualt values is 0.1
#   oncodrivefml_q (optional):   Q-value threshold for reporting OncodriveFML results. Defualt values is 0.001
#   oncodrivefml_conf (optional):   Directory and name of OncodriveFML configuration file
#   cgi (optional):   Name of folder and the results files from Cancer Genome Interpreter (CGI) analysis
#   clinical_info (optional):  Location of clinical data associated with each sample in MAF. Each file name (for each dataset) is expected to be separated by comma
#	  remove_duplicated_variants (optional):		Remove repeated variants in a particuar sample, mapped to multiple transcripts of same gene? Defulat value is "FALSE"
#   out_folder:   Name for the output folder that will be created within the directory with MAF files. If no output folder is specified the results will be saved in folder "Driver_analysis_report"
#   hide_code_btn: Hide the "Code" button allowing to show/hide code chunks in the final HTML report. Available options are: "TRUE" (default) and "FALSE"
#   ucsc_genome_assembly:  Human reference genome version used for signature analysis (default is "19")
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
  make_option("--maf_dir", action="store", default=NA, type='character',
              help="Directory with MAF files"),
  make_option("--maf_files", action="store", default=NA, type='character',
              help="List of MAF files to be processed"),
  make_option("--datasets", action="store", default=NA, type='character',
              help="Desired names of each dataset"),
  make_option("--samples_id_cols", action="store", default=NA, type='character',
              help="The name(s) of MAF file(s) column containing samples' IDs"),
  make_option("--dnds_p", action="store", default=0.05, type='double',
              help="dNdS method p-value threshold for reporting significant genes"),
  make_option("--dnds_q", action="store", default=1, type='double',
              help="dNdS method q-value threshold for reporting significant genes"),
  make_option("--activedriverwgs_p", action="store", default=0.05, type='double',
              help="ActiveDriverWGS method p-value threshold for reporting significant genes"),
  make_option("--activedriverwgs_fdr", action="store", default=1, type='double',
              help="ActiveDriverWGS method FDR value threshold for reporting significant genes"),
  make_option("--activedriverwgs_cores", action="store", default=1, type='double',
              help="Number of cores to be used for running ActiveDriverWGS method"),
  make_option("--activedriverwgs_all_genes", action="store", default=FALSE, type='logical',
              help="Run ActiveDriverWGS method for all genes"),
  make_option("--oncodriveclust_fdr", action="store", default=0.5, type='double',
              help="OncodriveClust method false discovery rate (FDR) threshold for reporting significant genes"),
  make_option("--ratios_ci", action="store", default=FALSE, type='logical',
              help="Calculate per-gene confidence intervals for the dN/dS ratios"),
  make_option("--hypermut_sample_cutoff", action="store", default=3000, type='integer',
              help="Mutations per gene to define ultra-hypermutator samples"),
  make_option("--max_muts_per_gene", action="store", default=3, type='integer',
              help="Maximum mutations per gene in same sample"),
  make_option("--out_folder", action="store", default="Driver_analysis_report", type='character',
              help="Output directory"),
  make_option("--genes_list", action="store", default="none", type='character',
              help="Location and name of a file listing genes of interest to be considered in the report"),
  make_option("--genes_blacklist", action="store", default="none", type='character',
              help="Location and name of a file listing genes to be excluded"),
  make_option("--samples_blacklist", action="store", default="none", type='character',
              help="Location and name of a file listing samples to be excluded"),
  make_option("--nonSyn_list", action="store", default=NA, type='character',
              help="List of variant classifications to be considered as non-synonymous"),
  make_option("--oncodrivefml", action="store", default="none", type='character',
              help="Name of folder and the results files from OncodriveFML analysis"),
  make_option("--oncodrivefml_p", action="store", default=0.01, type='double',
              help="P-value threshold for reporting OncodriveFML results. Defualt values is 0.1"),
  make_option("--oncodrivefml_q", action="store", default=0.1, type='double',
              help="Q-value threshold for reporting OncodriveFML results. Defualt values is 0.01"),
  make_option("--oncodrivefml_conf", action="store", default="none", type='character',
              help="Directory and name of OncodriveFML configuration file"),
  make_option("--cgi", action="store", default="none", type='character',
              help="Name of folder and the results files from CGI analysis"),
  make_option("--clinical_info", action="store", default="none", type='character',
              help="Location of clinical data associated with each sample in MAF"),
  make_option("--remove_duplicated_variants", action="store", default=TRUE, type='logical',
              help="Remove repeated variants in a particuar sample, mapped to multiple transcripts of same gene?"),
  make_option("--hide_code_btn", action="store", default=TRUE, type='logical',
              help="Hide the \"Code\" button allowing to show/hide code chunks in the final HTML report"),
  make_option("--ucsc_genome_assembly", action="store", default=19, type='integer',
              help="human reference genome version used for signature analysis")
)

opt <- parse_args(OptionParser(option_list=option_list))

##### Collect MAF files and correspondiong datasets names
opt$maf_files <- gsub("\\s","", opt$maf_files)
opt$clinical_info <- gsub("\\s","", opt$clinical_info)
opt$datasets <- gsub("\\s","", opt$datasets)

##### Read in argument from command line and check if all were provide by the user
if (is.na(opt$maf_dir) || is.na(opt$maf_files) || is.na(opt$datasets) ) {
  
  cat("\nPlease type in required arguments!\n\n")
  cat("\ncommand example:\n\nRscript driverAnalysis.R --maf_dir /data --maf_filessimple_somatic_mutation.open.PACA-AU.maf,PACA-CA.icgc.simple_somatic_mutation.maf --datasets ICGC-PACA-AU,ICGC-PACA-CA --dnds_p 0.05 --ratios_ci FALSE --hypermut_sample_cutoff 1000 --max_muts_per_gene 3 --ucsc_genome_assembly 19 --out_folder Driver_analysis_report\n\n")
  q()
  
} else if ( length(unlist(strsplit(opt$maf_files, split=',', fixed=TRUE))) != length(unlist(strsplit(opt$datasets, split=',', fixed=TRUE))) ) {
  
  cat("\nMake sure that the number of datasets names match the number of queried MAF files\n\n")
  q()
}

if ( !is.na(opt$samples_id_cols) && length(unlist(strsplit(opt$maf_files, split=',', fixed=TRUE))) != length(unlist(strsplit(opt$samples_id_cols, split=',', fixed=TRUE))) ) {
  
  cat("\nMake sure that the number of samples' ID columns match the number of queried MAF files\n\n")
  
  q()
}

##### Pre-define list of variant classifications to be considered as non-synonymous. Rest will be considered as silent variants. Default uses Variant Classifications with High/Moderate variant consequences (http://asia.ensembl.org/Help/Glossary?id=535)
if ( is.na(opt$nonSyn_list) ) {
  opt$nonSyn_list<- c("Frame_Shift_Del","Frame_Shift_Ins","Splice_Site","Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation")
} else {
  opt$nonSyn_list <- unlist(strsplit(opt$nonSyn_list, split=',', fixed=TRUE))
}

##### Set defualt paramters
if ( opt$ucsc_genome_assembly !=19 && opt$ucsc_genome_assembly !=38   ) {
  cat("\nCurrently human reference genome versions \"19\" and \"38\" are supported.\n\n")
  q()
}

if ( opt$ucsc_genome_assembly == 19 ) {
  ensembl_version <- 75
  
} else if ( opt$ucsc_genome_assembly == 38 ) {
  ensembl_version <- 86
}


##### Collect parameters

param_list <- list(maf_dir = opt$maf_dir,
                   maf_files = opt$maf_files,
                   datasets = opt$datasets,
                   samples_id_cols = opt$samples_id_cols,
                   dnds_p = opt$dnds_p,
                   dnds_q = opt$dnds_q,
                   activedriverwgs_p = opt$activedriverwgs_p,
                   activedriverwgs_fdr = opt$activedriverwgs_fdr,
                   activedriverwgs_cores = opt$activedriverwgs_cores,
                   activedriverwgs_all_genes = opt$activedriverwgs_all_genes,
                   oncodriveclust_fdr = opt$oncodriveclust_fdr,
                   ratios_ci = opt$ratios_ci,
                   hypermut_sample_cutoff = opt$hypermut_sample_cutoff,
                   max_muts_per_gene = opt$max_muts_per_gene,
                   ucsc_genome_assembly = opt$ucsc_genome_assembly,
                   out_folder = opt$out_folder,
                   genes_list = opt$genes_list,
                   genes_blacklist = opt$genes_blacklist,
                   samples_blacklist = opt$samples_blacklist,
                   nonSyn_list = opt$nonSyn_list,
                   oncodrivefml = opt$oncodrivefml,
                   oncodrivefml_p = opt$oncodrivefml_p,
                   oncodrivefml_q = opt$oncodrivefml_q,
                   oncodrivefml_conf = opt$oncodrivefml_conf,
                   cgi = opt$cgi, clinical_info = opt$clinical_info,
                   remove_duplicated_variants = opt$remove_duplicated_variants,
                   hide_code_btn = opt$hide_code_btn,
                   ucsc_genome_assembly = as.numeric(opt$ucsc_genome_assembly),
                   ensembl_version = as.numeric(ensembl_version)
)

##### Pass the user-defined argumentas to the summariseMAFs.R markdown script and run the analysis
rmarkdown::render(input = "driverAnalysis.Rmd",
                  output_dir = paste(opt$maf_dir, opt$out_folder, sep = "/"),
                  output_file = paste0(opt$out_folder, ".html"),
                  params = param_list)

##### Remove the assocaited MD file and the redundant folder with plots that are imbedded in the HTML report
unlink(paste0(opt$maf_dir, "/", opt$out_folder, "_files"), recursive = TRUE)

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
