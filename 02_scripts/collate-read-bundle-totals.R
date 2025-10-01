library(dplyr)
library(stringr)

files <- list.files("02_WG-RB-totals/", ".RB-totals-whole-genome.tsv", full.names = T)

rb_metrics_list <- list()

for (i in 1:length(files)) {
        rb_metrics <- read.table(files[i], header = T, row.names = 1)
        colnames(rb_metrics) <- str_replace_all(colnames(rb_metrics), "\\.", "-")
        rb_metrics_list[[i]] <- rb_metrics
}

rb_metrics_df <- dplyr::bind_cols(rb_metrics_list)

write.table(rb_metrics_df, "02_WG-RB-totals/read-bundle-totals-whole-genome.tsv", sep = "\t", quote = F, row.names = T, col.names = NA)
