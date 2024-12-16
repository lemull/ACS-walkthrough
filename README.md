Project Description

This repository contains data, code, and results for the study investigating the effects of Antenatal Corticosteroids (ACS) exposure on newborns' DNA methylation and their association with key phenotypes. The project involves preprocessing DNA methylation data, running Epigenome-Wide Association Studies (EWAS), and performing pathway and gene enrichment analyses to uncover significant biological patterns.

Data Preprocessing: Using Champ pipeline for raw data processing and probe-level extraction, and PCG pipeline for EWAS analysis.
Statistical Modeling: Linear mixed-effects models to evaluate the association between CpG methylation and phenotypes.
Pathway Analysis: Gene Ontology (GO) enrichment analysis and GSEA to identify biological pathways.
Phenotype Tables: Summarizing phenotypic data using descriptive statistics.
Time-point Comparisons: Methylation analysis at 1 month and 18 months after birth.

ACS/
│
├── Data/                      # Contains all initial raw data (metadata)
│
├── Code/  
│   ├── 00Select_Sample.R      # Select desir samples and filters samples for analysis: preterm (<37wk); received at least 1 treatment; have both visit (before 1 moth and after 18 month) 
│   ├── 01Sample_sheet.R       # Generates sample sheets for methylation analysis
│   ├── 02Copy_idat.R          # Copies and organizes idat files for selected samples
│   ├── 03modify_SampleSheet.R # Modifies and updates sample sheets (e.g. add covariate, metadata fixes)
│   ├── 04ChAMP.R              # Runs the ChAMP pipeline for preprocessing methylation data
│   ├── Descriptive_table.R    # Generates descriptive phenotype tables for the selected data
│   ├── Mix_effect_model.R     # Applies linear mixed-effects models to assess CpG-phenotype associations
│   └── PCG/                   # Contains PCG pipeline scripts for EWAS and GO analysis
│ 
├── Methylation_1/  
│   ├── Methylation_1.csv      # Sample sheet for Methylation_1 (1 month post-birth)  
│   └── idat_files/            # DNA methylation idat files for selected samples  
│
├── Methylation_18/  
│   ├── Methylation_18.csv     # Sample sheet for Methylation_18 (18 months post-birth)  
│   └── idat_files/            # DNA methylation idat files for selected samples  
│
├── Result/  
│   ├── ew_results/            # EWAS results (significant CpGs, beta values)  
│   ├── go_analysis/           # GO enrichment analysis outputs  
│   ├── gsea_results/          # GSEA results for biological pathways  
│   └── visualization/         # Plots and summary datasets  
│
└── README.md                  # Project documentation  
