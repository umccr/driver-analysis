# Driver analysis:

Ref: https://www.nature.com/articles/s41586-020-1965-x

Oncodrive has been superseded by OncodriveCLUSTL. See http://bg.upf.edu/group/projects/oncodrive-clust.php
Estimating background scores from synonymous variants. Look into updating code in the Rmd oncodriveCLUST_analysis chunk.

**dndscv:**

Ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5720395/

Reference coding sequence file CDS file is used as reference

R output
		
	[1] Loading the environment...
	
The p value for dndscv method is critical as it removes a lot of genes after filtering on this line:

	signif_genes_table.datasets[[i]]$dNdScv <- signif_genes_table.datasets[[i]]$dNdScv[ signif_genes_table.datasets[[i]]$dNdScv$pglobal_cv < params$dnds_p, ]
	

**ActiveDriverWGS:**

Larger WGS datasets of specific cancer types coupled with transcriptomes are required to find infrequent mutations in the vast non-coding genome and to outline their gene-regulatory functions. 
	
Quick tutorial http://htmlpreview.github.io/?https://github.com/reimandlab/ActiveDriverWGSR/blob/master/doc/ActiveDriverWGSR.html

Currently running analysis on cancer gene list only

			§ if ( !params$activedriverwgs_all_genes )

Consider running with `params$activedriverwgs_all_genes` set to true.

The current analysis uses ensembl db version 86 to grab genomic coordinates for key elements. Consider updating genomic version and using elements for non-coding regions
	
It takes a while to run the tool, even with the subset of data for cancers only genes.
	
The paper has ideas around using multiple methods for predicting drivers. Also look through supplementary material.

**OncoDriveCLUST:**

Oncodrive has been superseeded by OncodriveCLUSTL

		- https://bitbucket.org/bbglab/oncodriveclustl/src/master/
		- http://bg.upf.edu/group/projects/oncodrive-clust.php

		
Propose to construct the background model using the degree of clustering of synonymous mutations, which are assumed not to be under positive selection and may thus reflect the baseline mutation clustering of the tumour.

Detailed discussion on the installtion and usage of the tool can be found in [OncodriveCLUSTL.md](`./oncodriveCLUSTL/OncodriveCLUSTL.md`).

