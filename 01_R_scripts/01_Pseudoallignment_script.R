rm(list = ls())

# Working directories
raw_data_dir <- "/Volumes/AGAMBIAE_25/01.RawData"
output_dir <- "02_kallisto_output"
reference_genome <- "AnGambNW_F1_1_cds.idx"  

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# List all subdirectories in the raw data directory
subdirs <- list.dirs(raw_data_dir, full.names = TRUE, recursive = FALSE)

# Run kallisto for each subdirectory (sample)
for (subdir in subdirs) {
  sample <- basename(subdir)
  sample_files <- list.files(subdir, pattern = "*.fq.gz", full.names = TRUE)
  
  # Create a unique output directory for each sample
  sample_output_dir <- file.path(output_dir, sample)
  dir.create(sample_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Check if the number of files is even
  if (length(sample_files) %% 2 == 0) {
    # Paired-end reads
    kallisto_cmd <- paste(
      "kallisto quant",
      "-i", reference_genome,
      "-o", sample_output_dir,
      paste(sample_files, collapse = " ")
    )
  } else {
    # Single-end reads
    fragment_length_mean <- 150  # Actual fragment length mean
    fragment_length_sd <- 5      # Length standard deviation
    kallisto_cmd <- paste(
      "kallisto quant",
      "--single",
      "-i", reference_genome,
      "-o", sample_output_dir,
      "-l", fragment_length_mean,
      "-s", fragment_length_sd,
      paste(sample_files, collapse = " ")
    )
  }
  
  # Run kallisto command
  result <- system(kallisto_cmd, intern = TRUE)
  print(result)
}
