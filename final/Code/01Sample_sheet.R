# 加载必要的库
library(dplyr)
library(stringr)

# Step 1: 读取两个 CSV 文件
matched_df_1 <- read.csv("~/Desktop/matched_df_1.csv", stringsAsFactors = FALSE)
matched_df_18 <- read.csv("~/Desktop/matched_df_18.csv", stringsAsFactors = FALSE)

# Step 2: 加载行为数据并提取 study_corticosteroids 信息
load("~/Desktop/Data/Behavior.RData")

# 提取 src_subject_id 和 study_corticosteroids，用于后续判断 Sample_Group
study_info <- matched_df_18 %>%
  select(src_subject_id, study_corticosteroids, sex)

# Step 3: 定义函数来处理每个数据框
process_matched_data <- function(df, study_info) {
  df %>%
    left_join(study_info %>% select(src_subject_id, study_corticosteroids, sex), by = "src_subject_id") %>%
    # 合并 study_info 来获取 study_corticosteroids 信息
    left_join(study_info, by = "src_subject_id") %>%
    # 添加 Sample_Name
    mutate(Sample_Name = str_remove(patient_id_biorepository, "sub-"),
           # 根据 study_corticosteroids 确定 Sample_Group
           Sample_Group = ifelse(study_corticosteroids %in% c(1, 2), "Yes", "No"),
           # 固定 preterm 和 Sample_Plate
           preterm = 1,
           Sex = sex,
           Sample_Plate = NA,
           # extract Sentrix_ID and Sentrix_Position
           Sentrix_ID = str_extract(data_file1, "(?<=/)[0-9]{12}(?=_)"),
           Sentrix_Position = str_extract(data_file1, "(?<=_)[A-Z][0-9]{2}[A-Z][0-9]{2}(?=_|\\.idat)")
    ) %>%
    # 选择所需的列
    select(Sample_Name, Sample_Group, preterm, Sex, Sentrix_ID, Sentrix_Position, Sample_Plate, interview_age)
}


# 重新运行代码
processed_df_1 <- process_matched_data(matched_df_1, study_info)
write.csv(processed_df_1, "~/Desktop/1/Methylation_1/Methylation_1.csv", row.names = FALSE)
processed_df_18 <- process_matched_data(matched_df_18, study_info)
write.csv(processed_df_18, "~/Desktop/1/Methylation_18/Methylation_18.csv", row.names = FALSE)

