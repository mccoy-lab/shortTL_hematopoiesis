#!/usr/bin/env Rscript

# Plot chromosome-stratified depth distributions from mosdepth output

# Load libraries
library(tidyverse)
library(cowplot)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <-args[2]

df <- read_table(file=input_file) %>% 
    mutate(chrom = factor(chrom,levels=c("chr1","chr2","chr3","chr4","chr5",
                                       "chr6","chr7","chr8","chr9","chr10",
                                       "chr11","chr12","chr13","chr14","chr15",
                                       "chr16","chr17","chr18","chr19","chr20",
                                       "chr21","chr22")),
            sample = as.factor(sample))

g<-ggplot(data = df,
            aes(x=depth,y = proportion)) +
            geom_vline(xintercept=15, linetype="dashed", color = "royalblue") +
            geom_hline(yintercept=0.5, linetype="dashed", color = "royalblue") +
            geom_line(aes(group=sample),alpha=0.25) +
            scale_x_continuous(breaks=scales::pretty_breaks()) +
            facet_wrap(~chrom) +  
            theme_cowplot() +
            theme(strip.background = element_blank(),
                    legend.position="none") +
  labs(x="Depth",y="Proportion of sites")

ggsave(plot = g, file = output_file, width = 15, height = 10, dpi = 300)