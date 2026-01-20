#!/usr/bin/env Rscript

# Plot per-sample read counts from aggregated fastp data

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Read the data from command line argument
args <- commandArgs(trailingOnly = TRUE)
data <- read.csv(args[1], sep="\t")
output <- args[2]

# Extract Sample and Lane information from the file paths
data <- data %>%
    mutate(Sample = sub(".*/(.*_S\\d+)_.*", "\\1", filename),
           Lane = sub(".*_L(\\d+)_.*", "L\\1", filename),
           Reads = total_reads)

# Summarize total reads per sample
data_summary <- data %>%
    group_by(Sample) %>%
    summarize(Total_Reads = sum(Reads)) %>%
    ungroup() %>%
    arrange(Total_Reads)

# Merge the summarized data back to the original data
data <- data %>%
    inner_join(data_summary, by = "Sample")

# Create the stacked barplot
g<-ggplot(data, aes(x = reorder(Sample, Total_Reads), y = Reads/1e9)) +
    geom_bar(stat = "identity",color="white") +
    labs(x = "Sample",
         y = "Total Reads (billions)",
         fill = "Seq. run") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())

ggsave(g,filename = output, width = 15, height = 10, dpi = 300)