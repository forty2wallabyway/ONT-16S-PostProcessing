abund_files <- snakemake@input[["abund"]]
stats_files <- snakemake@input[["stats"]]

abund_out <- snakemake@output[["abund"]]
stats_out <- snakemake@output[["stats"]]

library(dplyr)

# Extract run names from file paths
get_run <- function(path) {
  basename(dirname(path))
}

# ---------------------------
# Merge abundance tables with run prefixes
# ---------------------------
abund_list <- lapply(abund_files, function(f) {
  
  run <- get_run(f)
  df <- read.delim(f, row.names = 1)
  
  # Prefix sample names with run
  colnames(df) <- paste0(run, "_", colnames(df))
  
  return(df)
})

# Merge all
abund_merged <- Reduce(function(x, y) {
  merge(x, y, by="row.names", all=TRUE)
}, abund_list)

rownames(abund_merged) <- abund_merged$Row.names
abund_merged$Row.names <- NULL

# Fill missing taxa
abund_merged[is.na(abund_merged)] <- 0

write.table(abund_merged, abund_out, sep="\t", quote=FALSE)

# ---------------------------
# Merge alignment stats
# ---------------------------
stats_list <- lapply(stats_files, function(f) {
  
  run <- get_run(f)
  df <- read.delim(f)
  
  df$sample <- paste0(run, "_", df$sample)
  df$run <- run
  
  return(df)
})

stats_merged <- bind_rows(stats_list)

write.table(stats_merged, stats_out, sep="\t",
            row.names=FALSE, quote=FALSE)
