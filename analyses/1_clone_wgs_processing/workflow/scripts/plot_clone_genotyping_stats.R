#!/usr/bin/env Rscript

# Script will generate plots to describe genotyping stats and sensitivity

#  Load libraries
library(tidyverse)
library(cowplot)
library(plyr)

args <- commandArgs(trailingOnly = TRUE)
stats_file <- args[1]
sensitivity_file <- args[2]
sensitivity_out <- args[3]
output_image <- args[4]

# Read in the data
stats.df <- read_tsv(stats_file)
sensitivity.df <- read_tsv(sensitivity_file)

# Clean df of uninformative columns and column names; reformat sensitivity df
stats.df <- stats.df %>% 
            select(-c(1,2,7,8,11,12,13))
colnames(stats.df) <- c("sample","nRefHom","nNonRefHom","nHets",
                        "nIndels","average_depth","nMissing")
sensitivity.df <- sensitivity.df %>% 
                    select(-c(1,2,7,8,11,12,13))
colnames(sensitivity.df) <- c("sample","nRefHom","nNonRefHom","nHets",
                                "nIndels","average_depth","nMissing")

# Subset stats to remove individual germline libraries that are
# excluded from sensitivity.df
stats.df <- stats.df %>% filter(sample %in% sensitivity.df$sample)

total <- sensitivity.df[sensitivity.df$sample == "consolidated_germline",]$nHets
sensitivity.df <- sensitivity.df %>% 
    mutate(total_germ_het = total,
            sensitivity = nHets/total_germ_het) %>%
    select(sample,sensitivity,average_depth)
sensitivity.df %>% 
    select(sample,sensitivity) %>% 
    write.table(file = sensitivity_out,sep = "\t",quote = FALSE,row.names = FALSE)

# Plotting stats
p1<-ggplot(stats.df, aes(x=reorder(sample,-average_depth), y=average_depth)) +
    geom_bar(stat="identity", fill="gray30",color="white") +
    geom_hline(yintercept=15, linetype="dashed", color = "gray98") +
    scale_y_continuous(limits=c(0,round_any(stats.df %>% 
                                                filter(sample == "consolidated_germline") %>% 
                                                pull(average_depth),
                                            accuracy=10,
                                            f=ceiling)),
                        breaks=seq(0,round_any(stats.df %>% 
                                                filter(sample == "consolidated_germline") %>% 
                                                pull(average_depth),
                                            accuracy=10,
                                            f=ceiling),by=10)) +
    theme_cowplot() +
    theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank()) +
    labs(x="Sample", y="Mean variant depth")

# Plotting sensitivity
xorder <- stats.df %>% arrange(-average_depth) %>% pull(sample)
sensitivity.df <- sensitivity.df %>%
    mutate(sample = factor(sample, levels = xorder))

p2<-ggplot(sensitivity.df, aes(x=sample, y=sensitivity)) +
    geom_bar(stat="identity", fill="gray70",color="white") +
    geom_hline(yintercept=0.5, linetype="dashed", color = "black") +
    scale_y_continuous(labels = scales::percent,breaks=seq(0,1,0.2)) +
    theme_cowplot() +
    theme(axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.5, size = 8)) +
    labs(x="Sample", y="Sensitivity")

p<-plot_grid(p1,p2,axis="lb",align="v",ncol=1,rel_heights=c(1,1))
ggsave(plot = p, file = output_image, width = 15, height = 10, dpi = 300)
