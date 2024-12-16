
install.packages("cardx")
install.packages("gtsummary")
library(gtsummary)
library(dplyr)


##### 1 month #####
samplesheet_data_1 <- read.csv("/Users/bingbing/Desktop/Methylation_1/Methylation_1.csv")

# Create a summary table
descriptive_table_1 <- samplesheet_data_1 %>%
  select(Sample_Group, Sex, baby_gest_at_birth, 
         baby_babyhc, baby_birth_weight, 
         mother_age1, mother_education1, mother_bmi) %>%
  tbl_summary(
    by = Sample_Group, # Grouping variable
    missing = "ifany", # Show missing data
    statistic = list(
      all_continuous() ~ "{mean} ({sd})", # Mean (SD) for continuous variables
      all_categorical() ~ "{n} ({p}%)"   # Count (percentage) for categorical variables
    ),
    digits = all_continuous() ~ 2        # Round continuous variables to 2 decimals
  ) %>%
  add_p() %>% # Add p-values for group comparisons
  modify_header(label ~ "**Variable**") %>% # Rename column headers
  modify_caption("**Descriptive Characteristics of Study Participants**") %>% 
  bold_labels()

# Print the table
descriptive_table_1




##### 18 month #####
samplesheet_data_18 <- read.csv("/Users/bingbing/Desktop/Methylation_18/Methylation_18.csv")

# Create a summary table
descriptive_table_18 <- samplesheet_data_18 %>%
  select(Sample_Group, Sex, baby_gest_at_birth, 
         baby_babyhc, baby_birth_weight, 
         mother_age1, mother_education1, mother_bmi) %>%
  tbl_summary(
    by = Sample_Group, # Grouping variable
    missing = "ifany", # Show missing data
    statistic = list(
      all_continuous() ~ "{mean} ({sd})", # Mean (SD) for continuous variables
      all_categorical() ~ "{n} ({p}%)"   # Count (percentage) for categorical variables
    ),
    digits = all_continuous() ~ 2        # Round continuous variables to 2 decimals
  ) %>%
  add_p() %>% # Add p-values for group comparisons
  modify_header(label ~ "**Variable**") %>% # Rename column headers
  modify_caption("**Descriptive Characteristics of Study Participants**") %>% 
  bold_labels()

# Print the table
descriptive_table_18
