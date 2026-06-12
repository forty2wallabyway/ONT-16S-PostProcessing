#############################################
# Merge wf-16s outputs + build metadata
# (depth inferred from directory structure)
#############################################

library(dplyr)

abund_files <- snakemake@input[["abund"]]
stats_files <- snakemake@input[["stats"]]

abund_out <- snakemake@output[["abund"]]
stats_out <- snakemake@output[["stats"]]
meta_out  <- snakemake@output[["meta"]]

# Helper functions
get_run <- function(path) basename(dirname(path))
get_depth <- function(path) basename(dirname(dirname(path)))  # SNS or DNS

#############################################
# 1. Merge abundance tables
#############################################

abund_list <- lapply(abund_files, function(f) {

  run <- get_run(f)
  depth_folder <- get_depth(f)

  df <- read.delim(f, row.names = 1, check.names = FALSE)

  # Create globally unique sample IDs
  colnames(df) <- paste0(depth_folder, "_", run, "_", colnames(df))

  return(df)
})

abund_merged <- Reduce(function(x, y) {
  merge(x, y, by = "row.names", all = TRUE)
}, abund_list)

rownames(abund_merged) <- abund_merged$Row.names
abund_merged$Row.names <- NULL

# Fill missing taxa with 0
abund_merged[is.na(abund_merged)] <- 0

write.table(abund_merged, abund_out, sep = "\t", quote = FALSE)

#############################################
# 2. Merge alignment stats
#############################################

stats_list <- lapply(stats_files, function(f) {

  run <- get_run(f)
  depth_folder <- get_depth(f)

  df <- read.delim(f)

  df$sample_id <- paste0(depth_folder, "_", run, "_", df$sample)
  df$run <- run
  df$depth <- ifelse(depth_folder == "SNS", "shallow", "deep")

  return(df)
})

stats_merged <- bind_rows(stats_list)

write.table(stats_merged, stats_out,
            sep = "\t", row.names = FALSE, quote = FALSE)

#############################################
# 3. Build metadata
#############################################

meta <- stats_merged %>%
  select(sample_id, run, depth)

# Extract animal ID from sample_id
# Format: SNS_run1_QS57 → QS57
meta$animal_id <- sub("^(SNS|DNS)_[^_]+_", "", meta$sample_id)

write.csv(meta, meta_out, row.names = FALSE)

#############################################
# 4. Optional sanity check
#############################################

cat("Pairing check:\n")
print(table(meta$animal_id, meta$depth))

#############################################
# END SCRIPT
#############################################
