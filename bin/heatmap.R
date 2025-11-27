#!/usr/bin/env Rscript
library(tidyverse)
library(optparse)

option_list = list(
  make_option(c("--tsv"), type='character', default = NULL, help="Output from SUMMARIZING_CLUSTERS process", metavar='character'),
  make_option(c("--output"), type='character', default = NULL, help="Output name - without extention", metavar='character')
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# data = read_tsv("/nfs/users/nfs_b/be5/scratch/testing_nextflow/5main_results/jf_extracted100_nrLitMari2_T_extracted100_k21_LC150_minCov150_clustered_nr0_clustered_putativeSLs.tsv")
data = read_tsv(opt$tsv)
numeric_cols <- c("entropy", "avg_coverage", "supporting_kmers", "`#_transcripts_in_cluster`")

normalize_score = function(x) {
  (x -mean(x, na.rm = TRUE))
}

normalize_data = data %>% 
  mutate(across(where(is.numeric), normalize_score, .names = "normalized_{.col}"))

heatmap_data = normalize_data %>%
  select(SL_ID, starts_with("normalized_")) %>%
  pivot_longer(cols = starts_with("normalized_"), names_to = "Metric", values_to = "Normalized_Value") %>%
  mutate(
    Metric = str_replace(Metric, "normalized_", ""),
    Metric = str_replace(Metric, "`#_transcripts_in_cluster`", "#_transcripts_in_cluster"),
    Metric = factor(Metric, levels = c("entropy", "avg_coverage", "supporting_kmers", "#_transcripts_in_cluster")))

heatmap_plot = ggplot(heatmap_data, aes(x = SL_ID, y = Metric, fill = Normalized_Value)) +
                        geom_tile() + 
                        scale_fill_gradient2(
                          low = "blue",
                          mid = "white",
                          high = "red", 
                          midpoint = 0,
                          name = "Normalized value") +
                        theme_minimal() + 
                        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))

output = opt$output + "_heatmap.png"
ggsave(output, heatmap_plot, dpi = 600)                      
# heatmap_plot
