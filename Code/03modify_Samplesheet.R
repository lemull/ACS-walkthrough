rm(list = ls())

library(dplyr)
# Load Behavior.RData
load("/Users/bingbing/Desktop/Data/Behavior.RData")
D <- trial_data_final_filtered

# select variables baby_babyhc & baby_birth_weight
baby_data <- D %>%
  select(src_subject_id, mother_age1, mother_education1, mother_bmi,baby_babyhc, baby_birth_weight, baby_gest_at_birth)%>%
  mutate(src_subject_id = sub("^sub-", "", src_subject_id))



##### 1 month #####
# Read Methylation_FV_T_Pre_data.csv
sample_data_1 <- read.csv("/Users/bingbing/Desktop/1/Methylation_1/Methylation_1.csv", stringsAsFactors = FALSE)
sample_data_1$Sentrix_ID <- as.character(sample_data_1$Sentrix_ID)
sample_data_1$Sample_Name <- as.character(sample_data_1$Sample_Name)

# Rename Sample_Name as src_subject_id for combination
#sample_data_1 <- sample_data_1 %>%
#  rename(src_subject_id = `Sample_Name`)
colnames(sample_data_1)[colnames(sample_data_1) == "Sample_Name"] <- "src_subject_id"
str(sample_data_1)

merged_data <- merge(sample_data_1, baby_data, by = "src_subject_id", all.x = TRUE)

# change src_subject_id back to Sample_Name
# merged_data <- merged_data %>%
#  rename(Sample_Name = src_subject_id)
colnames(sample_data_1)[colnames(sample_data_1) == "src_subject_id"] <- "Sample_Name"

# save csv file
write.csv(merged_data, "/Users/bingbing/Desktop/1/Methylation_1/Methylation_1.csv", row.names = FALSE)


##### 18 month #####
# Read Methylation_FV_T_Pre_data.csv
sample_data_18 <- read.csv("/Users/bingbing/Desktop/1/Methylation_18/Methylation_18.csv", stringsAsFactors = FALSE)

# Rename sample_Name as src_subject_id for combination
#sample_data_18 <- sample_data_18 %>%
#  rename(src_subject_id = Sample_Name)
colnames(sample_data_18)[colnames(sample_data_18) == "Sample_Name"] <- "src_subject_id"

merged_data <- merge(sample_data_18, baby_data, by = "src_subject_id", all.x = TRUE)

# change src_subject_id back to Sample_Name
#merged_data <- merged_data %>%
#  rename(Sample_Name = src_subject_id)
colnames(sample_data_18)[colnames(sample_data_18) == "src_subject_id"] <- "Sample_Name"

# save csv file
write.csv(merged_data, "/Users/bingbing/Desktop/Methylation_18/Methylation_18.csv", row.names = FALSE)

