rm(list=ls())
gc()
setwd <- "Users/bingbing/Desktop/1/Methylation_1"
idat_dir <- "/Users/bingbing/Desktop/1/Methylation_1"

# Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# Install ChAMP
# BiocManager::install("ChAMP", force = TRUE)
# BiocManager::install(c("GO.db", "minfi", "ChAMPdata", "Illumina450ProbeVariants.db", "sva", "IlluminaHumanMethylation450kmanifest", "limma", "RPMM", "DNAcopy", "preprocessCore", "impute", "marray", "wateRmelon", "goseq", "plyr", "GenomicRanges", "qvalue", "isva", "doParallel", "bumphunter", "quadprog", "shiny", "shinythemes", "plotly", "RColorBrewer", "DMRcate", "dendextend", "IlluminaHumanMethylationEPICmanifest", "matrixStats", "missMethyl", "combinat"), force = TRUE)
library("ChAMP")
sample_sheet <-  read.csv("/Users/bingbing/Desktop/1/Methylation_1/Methylation_1.csv", stringsAsFactors = FALSE)
sample_sheet$Sentrix_ID <- as.character(sample_sheet$Sentrix_ID)
sample_sheet$Basename <- file.path(idat_dir, sample_sheet$Sentrix_ID, sample_sheet$Sentrix_Position)

# load data through Champ
myLoad <- champ.load(idat_dir,arraytype="EPIC")
myLoad$pd$Slide <- as.character(myLoad$pd$Slide)



# myLoad <- cham.load(testDir)
# Check distribution of CpGs
CpG.GUI(arraytype="EPIC")
save(Anno, file = "/Users/bingbing/Desktop/result/FV/Anno_1.RData")
# Quality Control
champ.QC()
# QC.GUI()
QC.GUI(beta = myLoad$beta,arraytype="EPIC")

# Normalization
myNorm <- champ.norm(arraytype="EPIC",cores=8)
QC.GUI()
myNorm <- champ.norm(method="SWAN",arraytype="EPIC",cores=8)
#QC.GUI(beta = myNorm,arraytype="EPIC")

# Confounding Variable Test
myCombat <- champ.runCombat(beta = myNorm, pd = myLoad$pd)
champ.SVD(beta = myCombat, pd = myLoad$pd)

# DNA Methylation Diff
myDMP <- champ.DMP(beta = myCombat, pheno = myLoad$pd$Sample_Group)
# myDMP <- champ.DMP(beta = myNorm, pheno = myLoad$pd$Sample_Group, compare.group = c("Yes", "No"))
DMP.GUI()

# DNA Methylation Region
myDMR <- champ.DMR(beta=myCombat,pheno=myLoad$pd$Sample_Group,method="Bumphunter")
DMR.GUI(DMR=myDMR)

# Gene Set Enriched Analysis
myGSEA <- champ.GSEA(beta=myCombat,DMP=myDMP[[1]], DMR=myDMR, arraytype="EPIC",adjPval=0.05, method="fisher")
head(myGSEA$DMP)
head(myGSEA$DMR)


###### Save data ######
#' Input: path to the idats, Phenotype file
#' Output: RGset(RGChannelSet), Mset(MethylSet), GMset(GemomicMethylSet) in a list
library(minfi)

RGset <- read.metharray.exp(idat_dir, recursive = TRUE, verbose = TRUE, force = TRUE)
GMset <- preprocessNoob(RGset)
M_values <- log2(myLoad$beta / (1 - myLoad$beta))
beta_values <- myLoad$beta  # Extract Beta-values

save(RGset, file = "/Users/bingbing/Desktop/result/FV/RGset.RData")
save(GMset, file = "/Users/bingbing/Desktop/result/FV/GMset.RData")
write.csv(M_values, file = "/Users/bingbing/Desktop/result/SV/Mset_18.csv", row.names = TRUE)
write.csv(beta_values, file = "/Users/bingbing/Desktop/result/SV/Beta_values_18.csv", row.names = TRUE)
save(M_values, beta_values, file = "/Users/bingbing/Desktop/result/SV/beta$Mset_18.RData")

save(Anno, file = "/Users/bingbing/Desktop/result/FV/Anno_1.RData")
pd <- myLoad$pd
save(pd, file = "/Users/bingbing/Desktop/result/SV/pd_18.RData")
save(myLoad, file = "/Users/bingbing/Desktop/result/SV/myLoad_18.RData")
save(myNorm, file = "/Users/bingbing/Desktop/result/SV/myNorm_18.RData")
save(myCombat, file = "/Users/bingbing/Desktop/result/SV/myCombat_18.RData")

save(myDMP, file = "/Users/bingbing/Desktop/result/SV/myDMP_18.RData")
save(myDMR, file = "/Users/bingbing/Desktop/result/SV/myDMR_18.RData")

DMP <- myGSEA$DMP
DMR <- myGSEA$DMR

write.csv(DMP, file = "/Users/bingbing/Desktop/result/SV/DMP_18.csv", row.names = TRUE)
write.csv(DMR, file = "/Users/bingbing/Desktop/result/SV/DMR_18.csv", row.names = TRUE)

save(DMP, DMR, file = "/Users/bingbing/Desktop/result/SV/GSEA_DMP&DMR_18.RData")
save(myGSEA, file = "/Users/bingbing/Desktop/result/SV/myGSEA_18.RData")






###### DMP ScotterPlot ########
library(ggplot2)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Select sample_group Dmp dataset
dmp_table <- myDMP$`Yes_to_No`
write.csv(dmp_table, file = "DMP_table.csv")

# Select each DMP INFO
dmp_data <- myDMP[[1]]
dmp_data$CpG <- rownames(dmp_data)

# Select 每个 CpG 位点的染色体信息和位置
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
dmp_data <- merge(dmp_data, anno[, c("chr", "pos", "Name")], by.x = "CpG", by.y = "Name")

# Creat Scotterplot with a standard line
p <- ggplot(dmp_data, aes(x = pos, y = deltaBeta * 100, color = deltaBeta * 100)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -5, linetype = "dashed", color = "black") +
  facet_wrap(~ chr, scales = "free_x") +
  labs(title = "Scatterplot of Differentially Methylated CpG Sites",
       x = "Chromosome Position",
       y = "% Methylation Difference") +
  theme_minimal()
p
dir.create("G:/project/image", showWarnings = FALSE, recursive = TRUE)

# save as png
ggsave(filename = "G:/project/image/scatterplot1.png", plot = p, width = 30, height = 15, dpi = 1000, limitsize = FALSE)
