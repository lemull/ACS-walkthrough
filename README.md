## Project Description

This repository contains data, code, and results for the study investigating the effects of Antenatal Corticosteroids (ACS) exposure on newborns' DNA methylation and their association with key phenotypes. The project involves preprocessing DNA methylation data, running Epigenome-Wide Association Studies (EWAS), and performing pathway and gene enrichment analyses to uncover significant biological patterns.

# Summary
-  Data Preprocessing: Generate sample sheets（Code/01Sample_sheet.R）.
-  Data Processing： Using Champ pipeline for raw data processing （Rscript Code/04ChAMP.R） and probe-level extraction, and PCG pipeline for EWAS analysis（Scripts in Code/PCG/）.
-  Statistical Modeling: Linear mixed-effects models to evaluate the association between CpG methylation and phenotypes（Rscript Code/Mix_effect_model.R）.
-  Pathway Analysis: Gene Ontology (GO) enrichment analysis and GSEA to identify biological pathways.
-  Phenotype Tables: Summarizing phenotypic data using descriptive statistics.
-  Time-point Comparisons: Methylation analysis at 1 month and 18 months after birth.


# Workflow：
ACS/
│
├── Data/  
│   └── raw/                   # All initial raw data files
│
├── Code/  
│   ├── 00Select_Sample.R      # Selects and filters samples for analysis
│   ├── 01Sample_sheet.R       # Generates sample sheets for methylation analysis
│   ├── 02Copy_idat.R          # Copies and organizes idat files for selected samples
│   ├── 03modify_SampleSheet.R # Modifies and updates sample sheets (e.g., metadata fixes)
│   ├── 04ChAMP.R              # Runs the ChAMP pipeline for preprocessing methylation data
│   ├── Descriptive_table.R    # Generates descriptive phenotype tables for the selected data
│   ├── Mix_effect_model.R     # Applies linear mixed-effects models to assess CpG-phenotype associations
│   └── PCG/                   # Contains PCG pipeline scripts for EWAS and GO analysis
│
├── Methylation_1/             # DNA methylation data and sample sheet for 1-month samples
│   ├── Methylation_1.csv      # Sample sheet for Methylation_1
│   └── idat_files/            # DNA methylation idat files for selected samples
│
├── Methylation_18/            # DNA methylation data and sample sheet for 18-month samples
│   ├── Methylation_18.csv     # Sample sheet for Methylation_18
│   └── idat_files/            # DNA methylation idat files for selected samples
│
├── Result/                    # Outputs: tables, plots, and analysis results
│
└── README.md                  # Project documentation file
