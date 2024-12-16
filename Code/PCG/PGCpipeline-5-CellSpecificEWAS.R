########################################################################
### PGC Pipeline - Script 5 - Cell Specific EWAS
#########################################################################
rm(list=ls())
gc()

# Load the code to install packages
source("/Users/bingbing/Desktop/EPIC_QC/install_needed_packages.R") # PATH TO SOURCE FILE, this is an example path, change accordingly

#' Install packages if not already installed
bioc_packages <- c("TOAST", "RefFreeEWAS", "EpiDISH",
                    "limma")
cran_packages <- c("feather", "tibble", "tools",
                   "nnls", "data.table")

install_bioconductor_pkgs(pkgs = bioc_packages)
install_cran_pkgs(pkgs = cran_packages)

#' Run test
check_installed(pkgs = c(bioc_packages, cran_packages))

#' Load all packages, if needed
lapply(c(bioc_packages, cran_packages), require, character.only = TRUE)

#' Load beta values
# we need CpGs as rownames, so add the
# column that contains CpGs as rownames in the data
# In our data "V1" is the column containing CpGs
beta <- fread("/Users/bingbing/Desktop/noob_qcd_crossReactiveProbesRemoved_combat_CP_wcovar_age_ptsd_sex.csv", data.table = F)
beta <- column_to_rownames(beta, var = "V1")

## Load phenotype file with MethylationID (combination of Sentrix ID and Sentrix Position) as the 1st column
pheno <- read.csv("/Users/bingbing/Desktop/Pheno_QC_1.csv", row.names = 1)

# Now all the samples in beta and phenotype should match
# if samples are not matching, we will get the samples that
# are in phenotype file only
if(!all(colnames(beta) %in% rownames(pheno)) | !all(rownames(pheno) %in% colnames(beta))){
  beta <- beta[, which(colnames(beta) %in% rownames(pheno))]
  dim(beta)
  pheno <- pheno[rownames(pheno)%in%colnames(beta),]
}else message("Samples in both are matching")

# Now check if the columns in beta and pheno are in same order
# if not we can order them using the following code
table(colnames(beta) == rownames(pheno))
if(!all(colnames(beta) == rownames(pheno))){
  beta <- beta[, order(colnames(beta))]
  pheno <- pheno[order(rownames(pheno)), ]
  stopifnot(all(colnames(beta) == rownames(pheno))) # check the order again
}else message("Data is already ordered")

## Define Variable Names
study <- "ACS" # name of the study, e.g. "GTP", "DNHS" etc.
Var <- "Sample_Group" # name of the ptsd variable, coded as: cases = 1 and controls = 0
sex <- "Sex" # name of the sex variable


# We need to make the columns that are categories as factors
# As an example, we have gender/PTSD as categories,
# so we will make them as factors
fact_cols <- c(Var,sex)
pheno[fact_cols] <- lapply(pheno[fact_cols], factor)

# Now we need pheno and cell proportions separate
# to use them in the model
pheno_final <- pheno[, c(sex,Var)]
head(pheno_final)

cell_prop <- pheno[, c("CD8T", "CD4T", "NK",
                       "Bcell", "Mono", "Neu")] # cell estimations, change names if needed

all(rownames(pheno_final) == rownames(cell_prop)) # checking order, all should be TRUE
all(colnames(beta) == rownames(pheno_final))

# Function to get the significant CpGs
# Only those CpGs with p < 0.05
get_significant <- function(results, cutoff){
  filtered <-  lapply(results, function(x){
    x <- rownames_to_column(x, var = "CpGs")
    x[which(x['fdr'] < cutoff), ]
  })
  return(filtered)
}

# 使用基础 R 进行转换
pheno_final$Sample_Group <- ifelse(pheno_final$Sample_Group == 1, "case", "control")

# Run TOAST csDM function
design <- makeDesign(pheno_final, Prop = cell_prop)
design$formula

# fit model
fitted_model <- fitModel(design, as.matrix(beta))

colnames(design$design_matrix)
res_acs <- csTest(fitted_model, coef = "Sample_Group" ,cell_type = "joint")
res_acs_sig <- get_significant(results = res_acs, cutoff = 0.05)

# Save the output in the directory
save(res_acs, file = paste("/Users/bingbing/Desktop/ACS_csDM_TOAST.Rdata")) # THIS IS THE FILE TO BE SHARED WITH PGC-PTSD EWAS GROUP

view(res_acs)

get_significant <- function(results, cutoff) {
  if (!is.data.frame(results)) {
    stop("Input results must be a data frame.")
  }
  results <- rownames_to_column(results, var = "CpGs")
  significant_results <- results[results$fdr < cutoff, ]
  
  return(significant_results)
}

res_acs_sig <- get_significant(results = res_acs, cutoff = 0.05)


print(res_acs_sig)

save(res_acs_sig, file = paste("/Users/bingbing/Desktop/ACS_csDM_TOAST_fdr.Rdata")) # THIS IS THE FILE TO BE SHARED WITH PGC-PTSD EWAS GROUP


# ---------------------- END ------------------------------

# 加载必要的包
library(ggplot2)
library(data.table)

# 读取数据文件
pheno_data <- fread("/Users/bingbing/Desktop/Pheno_QC_2.csv")

# 绘制散点图
ggplot(pheno_data, aes(x = Sample_Group, y = Neu)) +
  geom_point() +
  labs(title = "Scatter plot of Sample_Group vs Neu",
       x = "Sample Group",
       y = "Neu") +
  theme_minimal()

pheno_data$Index <- 1:nrow(pheno_data)  # 为每一行添加一个索引
data_grouped <- pheno_data[, .(Mean_Neu = mean(Neu)), by = .(Sample_Group, Index)]

# 绘制折线图
ggplot(data_grouped, aes(x = Index, y = Mean_Neu, color = factor(Sample_Group))) +
  geom_line() +
  labs(title = "Line plot of Neu by Sample_Group",
       x = "Number",
       y = "Neu proportion",
       color = "Sample Group") +
  theme_minimal()
# ---------------------- END ------------------------------

# 加载必要的包
library(limma)
library(data.table)

# 读取beta值数据和phenotype数据
beta <- fread("/Users/bingbing/Desktop/noob_qcd_crossReactiveProbesRemoved_combat_CP_wcovar_age_ptsd_sex.csv", data.table = F)
beta <- column_to_rownames(beta, var = "V1")

pheno <- fread("/Users/bingbing/Desktop/Pheno_QC_2.csv")

# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)

# 假设 fitted_model$coefficients 包含你需要的 CpGs, beta 值, p 值或 FDR
results <- as.data.frame(fitted_model$coefficients)

# 假设结果表中包含以下列：'CpGs', 'Beta', 'FDR'
colnames(results) <- c("CpGs", "Beta", "FDR")

# 筛选出显著的 CpGs
significant_cpgs <- get_significant(results, cutoff = 0.05)

# 将显著的 CpGs 转换为长格式，以便使用 ggplot2 绘图
significant_long <- significant_cpgs %>%
  bind_rows() %>%
  gather(key = "Cell_Type", value = "Beta_Value", CD8T, CD4T, NK, Bcell, Mono, Neu)

# 添加 Sample_Group 信息
significant_long <- significant_long %>%
  left_join(pheno_data, by = "CpGs")

# 绘制折线图，不同的线表示不同的 Sample_Group
ggplot(significant_long, aes(x = CpGs, y = Beta_Value, color = Sample_Group)) +
  geom_line(aes(group = interaction(Sample_Group, Cell_Type))) +
  facet_wrap(~Cell_Type, scales = "free_y") +
  labs(title = "Beta Values for Significant CpGs Across Different Cell Types",
       x = "CpGs",
       y = "Beta Value",
       color = "Sample Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))




# 这里假设你已经有了 'results' 数据框，包含 CpGs, beta 值，和 fdr 值等
# 例如 'results' 可以是 fitModel 的输出，需要你根据你的实际数据调整

# 计算显著的 CpGs
significant_cpgs <- get_significant(results, cutoff = 0.05)

# 将显著的 CpGs 转换为长格式，以便使用 ggplot2 绘图
significant_long <- significant_cpgs %>%
  bind_rows() %>%
  gather(key = "Cell_Type", value = "Beta_Value", CD8T, CD4T, NK, Bcell, Mono, Neu)

# 添加 Sample_Group 信息
significant_long <- significant_long %>%
  left_join(pheno_data, by = "CpGs")

# 绘制折线图，不同的线表示不同的 Sample_Group
ggplot(significant_long, aes(x = CpGs, y = Beta_Value, color = Sample_Group)) +
  geom_line(aes(group = interaction(Sample_Group, Cell_Type))) +
  facet_wrap(~Cell_Type, scales = "free_y") +
  labs(title = "Beta Values for Significant CpGs Across Different Cell Types",
       x = "CpGs",
       y = "Beta Value",
       color = "Sample Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
