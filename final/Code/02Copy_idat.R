# Load necessary libraries
library(dplyr)


# Copy idat files to the destination folder
idat_source_dir <- "/Users/bingbing/Desktop/Methylation"
idat_dest_dir <- "/Users/bingbing/Desktop/1/Methylation_1"

# Read the cleaned sample data
sample_data_cleaned <- read.csv("/Users/bingbing/Desktop/1/Methylation_1/Methylation_1.csv", stringsAsFactors = FALSE)

for (i in 1:nrow(sample_data_cleaned)) {
  sentrix_id <- sample_data_cleaned$Sentrix_ID[i]
  sentrix_position <- sample_data_cleaned$Sentrix_Position[i]
  
  # Ensure Sentrix_ID and Sentrix_Position are not empty
  if (is.na(sentrix_id) || is.na(sentrix_position) || sentrix_id == "" || sentrix_position == "") {
    print(paste("Missing Sentrix_ID or Sentrix_Position at row:", i))
    next
  }
  
  # Construct file names
  idat_red_file <- paste0(sentrix_id, "_", sentrix_position, "_Red.idat")
  idat_grn_file <- paste0(sentrix_id, "_", sentrix_position, "_Grn.idat")
  
  # Construct full paths
  red_file_path <- file.path(idat_source_dir, idat_red_file)
  grn_file_path <- file.path(idat_source_dir, idat_grn_file)
  
  # Debugging: Print the generated paths
  print(paste("Checking Red file:", red_file_path))
  print(paste("Checking Green file:", grn_file_path))
  
  # Check if the files exist
  if (file.exists(red_file_path)) {
    file.copy(red_file_path, file.path(idat_dest_dir, idat_red_file), overwrite = TRUE)
  } else {
    print(paste("Red file not found:", red_file_path))
  }
  
  if (file.exists(grn_file_path)) {
    file.copy(grn_file_path, file.path(idat_dest_dir, idat_grn_file), overwrite = TRUE)
  } else {
    print(paste("Green file not found:", grn_file_path))
  }
}
# Remove columns ending in '.y' and rename columns ending in '.x'
# cleaned_data <- sample_data_cleaned %>%
#  select(-ends_with(".y")) %>%
#  rename_with(~ sub("\\.x$", "", .), ends_with(".x"))

# Save the cleaned data to a new CSV file
# write.csv(cleaned_data, "/Users/bingbing/Desktop/both_sample_data_cleaned.csv", row.names = FALSE)
