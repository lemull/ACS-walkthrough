load("/Users/bingbing/Desktop/result/SV/Anno_18.RData")
load("/Users/bingbing/Desktop/result/SV/myLoad_18.RData")
load("/Users/bingbing/Desktop/result/SV/myNorm_18.RData")
load("/Users/bingbing/Desktop/result/SV/myGSEA_18.RData")
load("/Users/bingbing/Desktop/result/SV/myCombat_18.RData")
load("/Users/bingbing/Desktop/result/SV/myDMP_18.RData")


brain_related_genes <- read.csv("/Users/bingbing/Desktop/genes_matrix_csv/rows_metadata.csv")
# 从 GSEA 结果提取显著基因列表
sig_genes_dmp <- unlist(strsplit(myGSEA$DMP$Genes[1], " "))
sig_genes_dmr <- unlist(strsplit(myGSEA$DMR$Genes[1], " "))
brain_genes <- brain_related_genes$gene_symbol
  
# 查找与脑相关基因的交集
brain_genes_in_dmp <- intersect(sig_genes_dmp, brain_genes)
brain_genes_in_dmr <- intersect(sig_genes_dmr, brain_genes)



# 假设使用 clusterProfiler 包进行 GO 富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
# 构建显著基因向量
gene_list <- sig_genes_dmp  # 使用你的显著基因列表
ego <- enrichGO(gene = gene_list,
                OrgDb = org.Hs.eg.db,  # 使用人类基因数据库
                keyType = "SYMBOL",
                ont = "BP",  # 生物过程
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2)

# 可视化结果
dotplot(ego, showCategory = 10)

# 提取 OR 和 P 值
brain_related_gsea <- myGSEA$DMP[grep("brain|neuron", myGSEA$DMP$Gene_List, ignore.case = TRUE), ]
print(brain_related_gsea)


sample_to_cpg <- myNorm # For group 1
sample_to_cpg <- myCombat # For group 18

cpg_to_genes <- data.frame(
  CpG = rownames(myDMP[[1]]),                  # CpG 名称
  Gene = as.character(myDMP[[1]]$gene)        # 基因信息
)

# 分解基因字符串为独立行
cpg_to_genes <- cpg_to_genes %>%
  # separate_rows(Gene, sep = ";") %>%
  filter(Gene %in% brain_related_genes$gene_symbol)

save(cpg_to_genes, file = "/Users/bingbing/Desktop/result/SV/cpg_to_genes_18.RData")


brain_gene_counts <- apply(sample_to_cpg, 2, function(sample_beta) {
  # 筛选 Beta 值大于阈值的 CpG
  active_cpgs <- rownames(sample_to_cpg)[sample_beta > 0.2]
  # 获取与 active_cpgs 相关联的基因
  active_genes <- unique(cpg_to_genes$Gene[cpg_to_genes$CpG %in% active_cpgs])
  # 返回样本中脑相关基因的数量
  length(active_genes)
})

# 检查 brain_gene_counts 是否正确
print(brain_gene_counts)




# 安装和加载必要包
library(dplyr)
library(tidyr)
library(lme4)
library(ggplot2)

# 加载数据
# 假设已经加载了以下对象:
# myLoad: 包含 CpG 甲基化水平和样本信息
# myDMP: 差异甲基化分析结果
# brain_related_genes: 包含大脑相关基因信息
# Anno: 包含 CpG 位点的注释信息

# 将 Anno 的 Annotation 转为 data.frame
anno_df <- as.data.frame(Anno$Annotation)

# 筛选感兴趣的列 (CpG 名称和其他信息)
anno_genes <- anno_df %>%
  filter(!is.na(CpG)) %>%
  select(CpG, M.index, U.index)

# 从 myDMP 提取显著 CpG 和基因关系
dmp_cpg_to_genes <- data.frame(
  CpG = rownames(myDMP[[1]]),
  Gene = as.character(myDMP[[1]]$gene)
)

# 分解基因字符串为独立行，并筛选大脑相关基因
cpg_to_genes <- dmp_cpg_to_genes %>%
  separate_rows(Gene, sep = ";") %>%
  filter(Gene %in% brain_related_genes$gene_symbol)

# 保存 CpG 到基因的映射
save(cpg_to_genes, file = "/Users/bingbing/Desktop/result/SV/cpg_to_genes_18.RData")

# 提取与目标基因（BDNF, MAPT, DISC1）相关的 CpG 位点
target_genes <- c("BDNF", "MAPT", "DISC1")
target_cpgs <- cpg_to_genes %>%
  filter(Gene %in% target_genes) %>%
  pull(CpG)

# 打印目标 CpG 位点
print(target_cpgs)

library(tibble)
# 筛选 target_cpgs 中的 CpG 位点
myNorm_target <- myNorm[rownames(myNorm) %in% target_cpgs, ]

# 转换为数据框格式并添加 CpG 和 Sample_Name 信息
df_target <- as.data.frame(myNorm_target) %>%
  rownames_to_column(var = "CpG") %>%  # 添加 CpG 信息
  pivot_longer(cols = -CpG, names_to = "Sample_Name", values_to = "beta_value")  # 转换为长格式

# 查看结果
head(df_target)

myLoad_pd <- as.data.frame(myLoad$pd)

# 合并 df_target 和 myLoad$pd
beta_long_18 <- df_target %>%
  left_join(myLoad_pd, by = "Sample_Name")  # 按 Sample_Name 合并

# 检查合并结果
head(beta_long_18)

write.csv(beta_long_18, file = "/Users/bingbing/Desktop/result/SV/beta_long_18.csv", row.names = FALSE)





beta_long_1 <- read.csv("/Users/bingbing/Desktop/result/FV/beta_long_1.csv")
beta_long_2 <- read.csv("/Users/bingbing/Desktop/result/SV/beta_long_18.csv")

# 假设 beta_long_1 和 beta_long_2 是两组数据
# 确保两组数据有一致的列名
head(beta_long_1)
head(beta_long_2)

# 合并两组数据
beta_long_combined <- rbind(beta_long_1, beta_long_2)
head(beta_long_combined)



# 包括随机斜率
random_effect_model <- lmer(
  beta_value ~ Sample_Group * interview_age + 
    (1 + interview_age | CpG) + (1 | Sample_Name),
  data = beta_long_combined
)
summary(random_effect_model)

# 提取随机效应
cpg_effects <- ranef(random_effect_model)$CpG

# 检查随机效应加载值
head(cpg_effects)

# 计算两个因子的均值
factor_mean_1 <- mean(cpg_effects[["(Intercept)"]], na.rm = TRUE)  # 随机截距的均值
factor_mean_2 <- mean(cpg_effects[["interview_age"]], na.rm = TRUE)  # 随机斜率的均值

# 将因子均值加入数据框
beta_long_combined <- beta_long_combined %>%
  mutate(
    random_factor_1 = factor_mean_1, 
    random_factor_2 = factor_mean_2
  )

# 查看数据框
head(beta_long_combined)

final_model <- lmer(beta_value ~ Sample_Group * interview_age + 
                      random_factor_1 + random_factor_2 + 
                      (1 | CpG) + (1 | Sample_Name),
                    data = beta_long_combined
)
summary(final_model)


# 构建不包含 group 的模型
model_no_group <- lmer(
  beta_value ~ interview_age + 
    (1 + random_factor_1 | CpG) + (1 + random_factor_2 | Sample_Name),
  data = beta_long_combined
)
summary(model_no_group)

# 比较模型
anova(final_model, model_no_group)




library(ggplot2)
summary_data <- beta_long_combined %>%
  group_by(Sample_Group, interview_age) %>%
  summarise(mean_beta = mean(beta_value, na.rm = TRUE),
            sd_beta = sd(beta_value, na.rm = TRUE), .groups = "drop")


ggplot(summary_data, aes(x = interview_age, y = mean_beta, color = Sample_Group, group = Sample_Group)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, size = 1) +  # 线性回归线
  labs(title = "CpG-Specific Group-Time Interaction on DNA Methylation",
       x = "Interview Age", y = "DNA Methyaltion Level") +
  theme_minimal()




# 分组：添加新列 age_group
beta_long_combined <- beta_long_combined %>%
  mutate(age_group = ifelse(interview_age < 1, "<1", ">=18"))

# 提取感兴趣基因的 CpG 位点
target_genes <- c("BDNF", "MAPT", "DISC1")
target_cpgs <- cpg_to_genes %>%
  filter(Gene %in% target_genes) %>%
  pull(CpG)

# 筛选数据
target_methylation <- beta_long_combined %>%
  filter(CpG %in% target_cpgs)

# 按基因和分组计算甲基化水平
gene_methylation <- target_methylation %>%
  left_join(cpg_to_genes, by = "CpG") %>%
  group_by(Gene, age_group, Sample_Group) %>%
  summarise(
    mean_beta = mean(beta_value, na.rm = TRUE),
    sd_beta = sd(beta_value, na.rm = TRUE),
    .groups = "drop"
  )

# 查看结果
print(gene_methylation)

# 可视化
ggplot(gene_methylation, aes(x = Sample_Group, y = mean_beta, fill = Sample_Group)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  geom_errorbar(aes(ymin = mean_beta - sd_beta, ymax = mean_beta + sd_beta), 
                width = 0.2, position = position_dodge(0.9)) +
  facet_grid(age_group ~ Gene, scales = "free_y") +  # 分面：按 age_group 和 Gene
  labs(title = "Methylation Levels by Gene, Group, and Age Group",
       x = "Sample Group", y = "DNA Methylation Level") +
  scale_fill_manual(values = c("No" = "lightblue", "Yes" = "pink")) +
  theme_minimal()

