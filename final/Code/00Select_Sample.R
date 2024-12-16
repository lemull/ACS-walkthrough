library(dplyr)
load("~/Desktop/Data/Behavior.RData")

#  select preterm baby 
filtered_data <- trial_data_final_filtered %>%
  filter(baby_gest_at_birth < 37)
head(filtered_data)


###### match data ######
library(dplyr)

# Step 1: load data from genomics_file 
genomics_data <- read.delim("/Users/bingbing/Desktop/Data/genomics_sample03.txt", header = TRUE, stringsAsFactors = FALSE)

# Step 2: extract rows corresponding to src_subject_id of filtered_data from genomics_data
genomics_filtered <- genomics_data[genomics_data$data_file1_type == "idat" & 
                                     genomics_data$src_subject_id %in% filtered_data$src_subject_id, ]


# Step 3: select the first interview_age <= 1 for each src_subject_id
df_age_leq_1 <- genomics_filtered %>%
  filter(interview_age <= 1) %>%
  group_by(src_subject_id) %>%
  slice_min(interview_age) %>%
  ungroup()

df_age_leq_1 <- df_age_leq_1 %>%
  distinct(src_subject_id, .keep_all = TRUE)
write.csv(df_age_leq_1, "~/Desktop/df_age_leq_1.csv", row.names = FALSE)

df_age_leq_1 <- read.csv("~/Desktop/df_age_leq_1.csv") 

# merge study_corticosteroids from trial_data_final_filtered with final_df
final_df_1 <- df_age_leq_1 %>%
  left_join(trial_data_final_filtered %>% select(src_subject_id, study_corticosteroids), 
            by = "src_subject_id")

head(final_df_1)

# Step 4: select the first interview_age >= 18 for each src_subject_id
df_age_geq_18 <- genomics_filtered %>%
  filter(interview_age >= 18) %>%
  group_by(src_subject_id) %>%
  slice_min(interview_age) %>%
  ungroup()

write.csv(df_age_geq_18, "~/Desktop/df_age_geq_18.csv", row.names = FALSE)

#### mannully select the rows with interview_age >= 18 #####
df_age_leq_18 <- read.csv("/Users/bingbing/Desktop/df_age_geq_18.csv")
View(df_age_leq_18)

# merge study_corticosteroids from trial_data_final_filtered with final_df
final_df_18 <- df_age_leq_18 %>%
  left_join(trial_data_final_filtered %>% select(src_subject_id, study_corticosteroids), 
            by = "src_subject_id")

head(final_df_18)




# Step 5: retain rows from df_age_leq_1 that have the same src_subject_id in df_age_geq_18
matched_df_1 <- final_df_1 %>%
  filter(src_subject_id %in% final_df_18$src_subject_id)
write.csv(matched_df_1, "~/Desktop/matched_df_1.csv", row.names = FALSE)
head(matched_df_1)

# Step 6: 从 df_age_geq_18 中保留在 df_age_leq_1 中有相同 src_subject_id 的行
matched_df_18 <- final_df_18 %>%
  filter(src_subject_id %in% final_df_1$src_subject_id)
write.csv(matched_df_18, "~/Desktop/matched_df_18.csv", row.names = FALSE)
head(matched_df_18)



# Step 6: 将 matched_df_1 添加到 matched_df_18 中，并按 src_subject_id 分组
final_df <- bind_rows(matched_df_1, matched_df_18) %>%
  arrange(sample_id_original)
# 查看结果
head(final_df)



