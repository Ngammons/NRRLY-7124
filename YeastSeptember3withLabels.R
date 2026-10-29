# Set working directory
setwd("/Users/nicolegammons/Desktop/Yeast_NovoGene_August2024")

# Load necessary library
library(dplyr)

# List of input and output files
input_files <- c("X1.txt", "X2.txt", "X3.txt", "X4.txt", "X5.txt")
output_files <- c("X1_output_file.txt", "X2_output_file.txt", "X3_output_file.txt", "X4_output_file.txt", "X5_output_file.txt")

# Loop through each file, process, and write the output
for (i in seq_along(input_files)) {
  
  # Read the text file into a data frame
  data <- read.table(input_files[i], header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  # Add a new column that joins columns 1 and 2 with an underscore
  data$new_column <- paste(data[[1]], data[[2]], sep = "_")
  
  # Move the new column to the start
  data <- data %>% select(new_column, everything())
  
  # Write the modified data to a new text file without the header
  write.table(data, output_files[i], sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Load data from files into a list
input_files <- c("X1_output_file.txt", "X2_output_file.txt", "X3_output_file.txt", "X4_output_file.txt", "X5_output_file.txt")
input_data <- lapply(input_files, function(file) read.table(file, header = FALSE, sep = "\t"))

# Function to filter and summarize data
filter_and_summarize <- function(data, dataset_name) {
  data <- na.omit(data)
  good_quality <- subset(data, V7 >= 80)
  bad_quality <- subset(data, V7 < 80)
  
  summary_stats <- data.frame(
    Dataset = dataset_name,
    Total = nrow(data),
    Good = nrow(good_quality),
    Bad = nrow(bad_quality)
  )
  
  return(list(good_quality = good_quality, bad_quality = bad_quality, summary_stats = summary_stats))
}

# Apply filtering to all datasets and collect summary statistics
filtered_data <- lapply(seq_along(input_data), function(i) {
  filter_and_summarize(input_data[[i]], paste0("X", i))
})

# Combine all summary statistics into one data frame
summary_stats_df <- do.call(rbind, lapply(filtered_data, function(x) x$summary_stats))

# Save the summary statistics to a file
write.table(summary_stats_df, file = "Quality_Summary_Statistics.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Quality summary statistics saved to Quality_Summary_Statistics.txt\n")

# Extract good quality variants for X1
good_quality_X1 <- filtered_data[[1]]$good_quality

# Function to find unique variants not in X1
find_unique_variants <- function(data, reference) {
  return(data[!duplicated(rbind(data, reference))[(nrow(data)+1):nrow(rbind(data, reference))],])
}

# Identify unique variants in X2-X5 compared to X1
unique_variants_list <- lapply(2:5, function(i) {
  dataset <- filtered_data[[i]]$good_quality
  find_unique_variants(dataset, good_quality_X1)
})

# Find common V1 values across X2 to X5
common_V1_values <- Reduce(intersect, lapply(unique_variants_list, function(data) data$V1))

# Compile common variants across datasets
common_variants <- unique_variants_list[[1]][unique_variants_list[[1]]$V1 %in% common_V1_values, ]
for (i in 2:4) {
  common_data <- unique_variants_list[[i]][unique_variants_list[[i]]$V1 %in% common_V1_values, ]
  common_variants <- merge(common_variants, common_data, by = "V1", all = TRUE)
}

# Save common variants to file
write.table(common_variants, file = "Common_Good_Quality_Variants_Across_X2_to_X5.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("Common good quality variants saved to Common_Good_Quality_Variants_Across_X2_to_X5.txt\n")

# Find and save unique not common variants for each dataset with dataset origin column
all_not_common_variants <- do.call(rbind, lapply(2:5, function(i) {
  dataset <- unique_variants_list[[i-1]]
  unique_V1_values <- setdiff(dataset$V1, common_V1_values)
  unique_data <- dataset[dataset$V1 %in% unique_V1_values, ]
  unique_data <- na.omit(unique_data)
  
  # Add a column to identify the dataset origin
  unique_data$Dataset_Origin <- paste0("X", i)
  
  # Save to individual file
  output_file <- paste0("X", i, "_Unique_Good_Quality_Variants_Not_Common.txt")
  write.table(unique_data, file = output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat(paste("Unique good quality variants for X", i, "saved to", output_file, "\n"))
  
  return(unique_data)
}))

# Save consolidated unique variants across all datasets with the origin column
write.table(all_not_common_variants, file = "All_Not_Common_Good_Quality_Variants_Across_X2_to_X5.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("Consolidated not common good quality variants saved to All_Not_Common_Good_Quality_Variants_Across_X2_to_X5.txt with dataset origin column\n")

# Separate column 12 by | and expand it into multiple columns while retaining the Dataset_Origin column
separate_column_12 <- function(data) {
  if (ncol(data) >= 12) {
    split_data <- strsplit(as.character(data$V9), split = "\\|")
    max_len <- max(sapply(split_data, length))
    
    expanded_columns <- do.call(rbind, lapply(split_data, function(x) {
      length(x) <- max_len
      return(x)
    }))
    
    colnames(expanded_columns) <- paste0("V9_", seq_len(max_len))
    
    # Retain Dataset_Origin column and bind with expanded columns
    data <- cbind(data[, -12], expanded_columns, Dataset_Origin = data$Dataset_Origin)
  }
  return(data)
}

# Apply the separation to the consolidated not common variants
all_not_common_variants <- separate_column_12(all_not_common_variants)

# Save the modified table with separated column 12 and Dataset_Origin
write.table(all_not_common_variants, file = "All_Not_Common_Good_Quality_Variants_Across_X2_to_X5_Separated.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("Consolidated not common good quality variants with separated column 12 and Dataset_Origin saved to All_Not_Common_Good_Quality_Variants_Across_X2_to_X5_Separated.txt\n")

# Extract unique V1 values across all datasets
all_V1_values <- unique(unlist(lapply(filtered_data, function(x) x$good_quality$V1)))

# Count the total number of unique V1 values
total_unique_V1 <- length(all_V1_values)

# Save the total count to Unique_Locus_Tags.txt
write.table(total_unique_V1, file = "Unique_Locus_Tags.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("Total unique locus tags count saved to Unique_Locus_Tags.txt\n")
