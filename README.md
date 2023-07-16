# __GRN Reconstruction__
This repository contains code that was developed as part of the Bioinformatics & Biocomplexity MSc Thesis by Abdool Al-Khaledi at Utrecht University. The pipeline extends the capabilities of the DoRothEA network [1] and DecoupleR’s Multivariate Linear Model (MLM) [2] to generate cell type-specific Gene Regulatory Networks (GRNs). The generated networks can subsequently be employed for downstream analysis.

In this case, two separate analyses were utilized to investigate Autism Spectrum Disorder (ASD) related genes. The first method was a distributional analysis of engineered metrics followed by a Gene Set Enrichment Analysis (GSEA). The second method involved an Ordinary Least Squares (OLS) regression analysis to model the ASD network subset as a function of total network size and time.

# **Repository Contents**
Mouse orthologues of SFARI genes: This directory contains the "autism_genes.csv" file, which lists ASD-related genes curated by The Simons Foundation Autism Research Initiative (SFARI) [3], as well as the mouse orthologues of these genes.

GRN_reconstruction_pipeline: This script is responsible for the cell type-specific reconstruction of GRNs. The pipeline requires a clustered single cell/nucleus RNA-seq CellXGene matrix as input and outputs cell type-specific networks averaged across all time points (as network .csv files). If the age-delineated method is used, a folder for each timepoint is created containing the cell type-specific networks of the respective timepoint. The generated networks can be used in downstream analysis.

Network_distributional_analysis: This notebook utilizes whole data networks, although individual timepoint networks can also be used. The analysis quantifies metrics to measure GRN activity and ASD related network activity, identifies cell type-specific ASD-related enrichment, and finally, performs GSEA on the ASD TF-gene subsets for the enriched cohort.

Network_regression_analysis: This notebook makes use of age-delineated networks. The input should be a folder containing subdirectories for each timepoint, where each subdirectory contains cell-type specific GRNs. The notebook performs an OLS regression to delineate the effect of time on ASD-related GRN activity.

**Note: The notebooks can be downloaded if GitHub is unable to render them.**

**References**
1. Garcia-Alonso L, Holland CH, Ibrahim MM, Turei D, Saez-Rodriguez J. Benchmark and integration of resources for the estimation of human transcription factor activities. Genome Res. 2019 Aug;29(8):1363-1375. doi: 10.1101/gr.240663.118. Epub 2019 Jul 24. Erratum in: Genome Res. 2021 Apr;31(4):745. PMID: 31340985; PMCID: PMC6673718.
2. Badia-I-Mompel P, Vélez Santiago J, Braunger J, Geiss C, Dimitrov D, Müller-Dott S, Taus P, Dugourd A, Holland CH, Ramirez Flores RO, Saez-Rodriguez J. decoupleR: ensemble of computational methods to infer biological activities from omics data. Bioinform Adv. 2022 
3. Banerjee-Basu S, Packer A. SFARI Gene: an evolving database for the autism research community. Dis Model Mech. 2010 Mar-Apr;3(3-4):133-5. doi: 10.1242/dmm.005439. PMID: 20212079.
