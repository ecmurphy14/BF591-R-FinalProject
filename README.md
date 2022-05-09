# BF591-R-FinalProject

This repo hosts the code and csv files for Liz Murphy's R Shiny application final project for BF591-R.

This app explores data from GSE64810, accessible here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810.
All data files to run the app are available in this repository, in the folder final_project.
  The data folder includes:
    *sample_info.csv* processed from the SOFT files.
    
    *GSE64810_mlhd_DESeq2_norm_counts_adjust.tsv* normalized count matrix downloaded from the GEO.
    
    *GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.tsv* DESeq2 differential expression analysis results.
    
    *fgsea_results.csv* GSEA results as generated from the DESeq2 csv in processing.R
    

processing.R is an script with code for reading the sample metadata from GSE64810 SOFT files, as well as code for running fgsea on the available DESeq2 results.

app.R is the code for the application, available in the final_project folder. 
