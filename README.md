# __GRN Reconstruction__
This repository contains code that was developed as part of the Bioinformatics & Biocomplexity MSc Thesis by Abdool Al-Khaledi at Utrecht University. The pipeline extends the capabilities of the DoRothEA network and DecoupleR’s Multivariate Linear Model (MLM) to generate cell type-specific Gene Regulatory Networks (GRNs). The generated networks can subsequently be employed for downstream analysis.

In this case, two separate analyses were utilized to investigate Autism Spectrum Disorder (ASD) related genes. The first method was a distributional analysis of engineered metrics followed by a Gene Set Enrichment Analysis (GSEA). The second method involved an Ordinary Least Squares (OLS) regression analysis to model the ASD network subset as a function of total network size and time.

**Repository Contents**
Mouse orthologues of SFARI genes: This directory contains the "autism_genes.csv" file, which lists ASD-related genes curated by The Simons Foundation Autism Research Initiative (SFARI), as well as the mouse orthologues of these genes.

GRN_reconstruction_pipeline: This script is responsible for the cell type-specific reconstruction of GRNs. The pipeline requires a clustered single cell/nucleus RNA-seq CellXGene matrix as input and outputs cell type-specific networks averaged across all time points (as network .csv files). If the age delineated method is used, a folder for each timepoint is created containing the cell type-specific networks of the respective timepoint. The generated networks can be used in downstream analysis.

Network_distributional_analysis: This notebook utilizes whole data networks, although individual timepoint networks can also be used. The analysis quantifies metrics to measure GRN activity and ASD related network activity, identifies cell type-specific ASD related enrichment, and finally, performs GSEA on the ASD TF-gene subsets for the enriched cohort.

Network_regression_analysis: This notebook makes use of the age delineated networks. The input should be a folder containing subdirectories for each timepoint, where each subdirectory contains cell-type specific GRNs. The notebook performs an OLS regression to delineate the effect of time on ASD related GRN activity.
