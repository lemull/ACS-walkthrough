#########################################################################
### PGC Pipeline - Script 4 - EWAS
#########################################################################
#' Clean
rm(list=ls())
gc()

#' Load packages
load_pkgs <- c("CpGassoc", "data.table", "tibble", "feather")
lapply(load_pkgs, require, character.only = TRUE)

## Load Methylation Data
beta <- fread("/Users/bingbing/Desktop/noob_qcd_crossReactiveProbesRemoved_combat_CP_wcovar_age_acs_sex.csv", data.table = F)
rownames(beta)<-beta$V1
beta<-beta[,-1]

## Load phenotype file with MethylationID (combination of Sentrix ID and Sentrix Position) as the 1st column
pheno <- read.csv("/Users/bingbing/Desktop/Pheno_QC_1.csv", row.names = 1)

# A function here to get and order the required data
clean_order <- function(beta, pheno){
  cpg <- beta[, colnames(beta) %in% row.names(pheno)]
  cpg <- cpg[, order(colnames(cpg))]
  pheno <- pheno[rownames(pheno) %in% colnames(cpg), ]
  pheno <- pheno[order(rownames(pheno)), ]
  print(table(rownames(pheno) == colnames(cpg))) # should be TRUE
  stopifnot(all(rownames(pheno) == colnames(cpg)))
  return(list(pheno = pheno, cpg = cpg))
}

cleaned_df <- clean_order(beta = beta, pheno = pheno)
pheno <- cleaned_df$pheno
write.csv(pheno, file = "/Users/bingbing/Desktop/Pheno_QC_EWAS.csv",row.names = T)

## Define variables
study <- "ACS" # name of the study, e.g. "GTP", "DNHS" etc.
Var <- "Sample_Group" # name of the acs variable, coded as: cases = 1 and controls = 0
acs <- pheno[, Var, FALSE]

## acs variable should be numeric
cleaned_df$pheno[, Var] <- as.numeric(cleaned_df$pheno[, Var])

## Define covariates to be adjusted for EWAS
## Covariates to be included:
##  - cell types from step 3 ("CD8T","CD4T","NK","Bcell","Mono")
##  - GWAS PCs PC1 and PC2 (if available), if not mPC2 (Comp.2) and mPC3 (Comp.3) from step 3.1
##  - age
##  - sex (if applicable)
covar <- data.frame(pheno[,c("CD8T","CD4T","NK","Bcell","Mono","Neu","Sex")])

# Function to run EWAS with CpGAssoc
cpg_assoc_test <- function(cpg, pheno, covar){
  message("Running test, patience ...")
  test <- cpg.assoc(cpg, pheno, covar, logit.transform = T, large.data=TRUE)
  assoc <- test$results
  eff <- test$coefficients
  results <- cbind(assoc, eff)
  return(list(results = results, test = test))
}


# Run EWAS with CpGAssoc
results <- cpg_assoc_test(cpg = cleaned_df$cpg,
                          pheno = cleaned_df$pheno[, Var],
                          covar = covar)

save(results, file = paste0("/Users/bingbing/Desktop/", study,"_EWAS.Rdata")) # THIS IS THE FILE TO BE SHARED WITH PGC-acs EWAS GROUP

view(results)
