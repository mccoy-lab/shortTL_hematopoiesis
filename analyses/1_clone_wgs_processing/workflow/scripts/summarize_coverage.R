#!/usr/bin/env Rscript

# Generate a plot and tsv that summarized the per-sample mean coverage

#Load libraries
library(tidyverse)
library(cowplot)
library(ggbeeswarm)

#Read input
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
individual <- args[2]
output_tsv <- args[3]
output_plot <- args[4]

df <- read_table(input_file) %>% 
    mutate(consolidated_germline = ifelse(File == "consolidated_germline.bam.stats",TRUE,FALSE))
df$library <- strsplit(df$File,"\\.") %>% sapply(function(x) x[1])
df <- df %>% group_by(library)

#Compute mean coverage and save to tsv file
cov_mean.df <- lapply(unique(df$File), function(x) {
    df_sub <- df %>% filter(File == x)
    mean <- sum(df_sub$Cov*df_sub$Count)/sum(df_sub$Count)
    return(data.frame("library" = strsplit(x,"\\.")[[1]][1],
                      "mean_cov" = mean))
}) %>% 
    bind_rows()

cov_mean.df %>% write.table(file = output_tsv,sep = "\t",quote = FALSE,row.names = FALSE)

#Generate coverage plot and save to file
g1<-ggplot(data = df,aes(x=Cov,y=Count/1e6)) +
    geom_vline(xintercept = seq(10,100,by=10),
                color="gray95",
                linewidth=0.1) +
    geom_vline(xintercept = 15,
                color="gray50",
                linetype="dashed",
                linewidth=0.1) +
    geom_vline(data = cov_mean.df,
                aes(xintercept = mean_cov),
                color = "royalblue",
                linewidth=0.1) +
    geom_text(data = cov_mean.df, 
                aes(x = mean_cov + 3),
                y = (max(df$Count)/1e6)*0.99, 
                label = "bar(x)",
                parse=TRUE,
                color = "royalblue",
                vjust = 1, 
                hjust = 0, 
                size = 3) +
    geom_col(fill="black",
                color="black",
                linewidth=0.05) +
    facet_wrap(~library) +
    scale_x_continuous(limits = c(0,100),
                        breaks=seq(0,100,by=20)) +
    theme_cowplot() +
    theme(strip.text = element_text(size = 6),
            axis.text = element_text(size=6),
            axis.title = element_text(size=8,face="bold"),
            plot.title = element_text(size=10),
            strip.background = element_blank()) +
    labs(title = paste("Distribution of sequencing coverage per site: ",individual), 
            x = "Coverage/site", 
            y = "Count (millions)")

g2<-ggplot(data = cov_mean.df,
            aes(x="",y = mean_cov)) +
            geom_violin(fill = "royalblue",
                        color = "black",
                        alpha = 0.5,
                        linewidth=0.5) +
            geom_boxplot(fill="white",
                            width=0.3,
                            linewidth=0.5,
                            outlier.color="transparent") +
            geom_quasirandom(width=0.1,
                                size=0.5,
                                color="black",
                                alpha=1) +
            scale_y_continuous(limits = c(0,max(max(cov_mean.df$mean_cov),100)), 
                                breaks = seq(0,max(max(cov_mean.df$mean_cov),100),10),
                                position="right") +
            theme_cowplot() +
            theme(axis.text = element_text(size=6),
                    axis.text.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.line.x = element_blank(),
                    axis.title = element_text(size=8,face="bold"),
                    plot.title = element_text(size=10)) +
            labs(y = "Mean coverage/site")

g <- plot_grid(g1,g2,ncol=2,rel_widths = c(5,1))
ggsave(plot = g, file = output_plot, width = 15, height = 10, dpi = 300)