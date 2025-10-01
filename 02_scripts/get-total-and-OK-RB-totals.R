# get-total-and-OK-RB-totals.R

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript get-total-and-OK-RB-totals.R <input.RBs> <output.tsv>")
}

input_file <- args[1]
output_file <- args[2]

cat("Reading file:", input_file, "\n")
cat("Output file:", output_file, "\n")
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)


if (!file.exists(input_file)) {
  stop(paste("ERROR: Input file not found:", input_file))
}

library(stringr)
library(dplyr)

rbs_bck <- read.table(input_file, sep = "\t", header = FALSE, row.names = 1)
cat("Loaded RBs. Dimensions:", dim(rbs_bck), "\n")

if (ncol(rbs_bck) < 2) {
  stop("ERROR: RBs table must have at least 2 data columns.")
}

# Add capped values for visualization/stats
rbs_bck$x_tmp <- pmin(rbs_bck[, 1], 10)
rbs_bck$y_tmp <- pmin(rbs_bck[, 2], 10)
rbs_bck$x <- rbs_bck[,1]
rbs_bck$y <- rbs_bck[,2]

OK_RBS <- sum(rbs_bck$x >= 2 & rbs_bck$y >= 2)
TOTAL_RBS <- nrow(rbs_bck)

cat("OK_RBS:", OK_RBS, "TOTAL_RBS:", TOTAL_RBS, "\n")

RBs_df <- data.frame(
  sample = str_split(input_file, "\\.", simplify = TRUE)[,1],
  OK_RBS = OK_RBS,
  TOTAL_RBS = TOTAL_RBS
)

rownames(RBs_df) <- RBs_df$sample
RBs_df$sample <- NULL

RBs_df <- t(RBs_df)

write.table(RBs_df, output_file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
cat("Output written to:", output_file, "\n")
