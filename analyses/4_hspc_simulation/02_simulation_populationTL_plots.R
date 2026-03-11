# Plot simulation output

library(tidyverse)
library(data.table)
library(scales)
library(khroma)
library(patchwork)

#============#
# Data setup #
#============#

data <- fread("popTL_unif_simulations.txt") %>%
  .[, popsize_tertp := popsize * vaf] %>%
  .[is.na(popsize_tertp), popsize_tertp := 0]
data$mut <- factor(data$mut,
                   levels = unique(data$mut),
                   labels = c("1e-04", "1e-05", "1e-06", "1e-07", "1e-08", "1e-09"))

# Annotate with TL at birth
data2 <- data %>%
  setorder(mut, trial, year) %>%
  group_by(mut, trial) %>%
  mutate(birth_tl = wt_tl[1]) %>%
  setDT() %>%
  .[, ident := paste0(mut, "_", trial)]

data_summary2 <- data2[birth_tl >= 8500] %>%
  setDT() %>%
  .[, is_chip := vaf >= 0.01] %>%
  group_by(mut, year, id) %>%
  summarize(mean_vaf = mean(vaf),
            sd_vaf = sd(vaf),
            q1 = quantile(vaf, 0.05, na.rm = TRUE),
            q2 = quantile(vaf, 0.95, na.rm = TRUE),
            perc_chip = sum(is_chip) / n(),
            mean_tl_wt = mean(wt_tl, na.rm = TRUE),
            q1_tl_wt = quantile(wt_tl, 0.05, na.rm = TRUE),
            q2_tl_wt = quantile(wt_tl, 0.95, na.rm = TRUE),
            mean_tl_tertp = mean(tertp_tl, na.rm = TRUE),
            q1_tl_tertp = quantile(tertp_tl, 0.05, na.rm = TRUE),
            q2_tl_tertp = quantile(tertp_tl, 0.95, na.rm = TRUE),
            mean_popsize_wt = mean(popsize),
            q1_popsize_wt = quantile(popsize, 0.05, na.rm = TRUE),
            q2_popsize_wt = quantile(popsize, 0.95, na.rm = TRUE),
            mean_popsize_tertp = mean(popsize_tertp),
            q1_popsize_tertp = quantile(popsize_tertp, 0.05, na.rm = TRUE),
            q2_popsize_tertp = quantile(popsize_tertp, 0.95, na.rm = TRUE),
            mean_lineages = mean(n_tert_lineages),
            q1_lineages = quantile(n_tert_lineages, 0.05, na.rm = TRUE),
            q2_lineages = quantile(n_tert_lineages, 0.95, na.rm = TRUE)) %>%
  setDT()

fwrite(data2, "popTL_sim_data.txt", sep = "\t")


#=============================#
# Plots of all mutation rates #
#=============================#

# subset to mutation rates of interest
data_plot <- data2 %>%
  .[mut %in% c("1e-04", "1e-06", "1e-08")]
data_plot$mut <- factor(data_plot$mut,
                        levels = c("1e-04", "1e-06", "1e-08"),
                        labels = c("italic(mu) == 1 %*% 10^-4",
                                   "italic(mu) == 1 %*% 10^-6",
                                   "italic(mu) == 1 %*% 10^-8"))

p1 <- ggplot(data_plot[mut == "italic(mu) == 1 %*% 10^-4"],
       aes(x = year, y = vaf)) +
  geom_line(aes(group = trial, color = birth_tl), linewidth = 0.3) +
  theme_bw() +
  xlim(0, 90) +
  ylim(0, 0.0413) +
  facet_wrap(~ mut, scales = "free", labeller = label_parsed) +
  scale_color_viridis_c(option = "plasma") +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "gray90", color = "black"),
        strip.text = element_text(size = 10, margin = margin(t = 2, b = 2)),
        text = element_text(family = "Helvetica"),
        legend.position = "none") +
  labs(x = "Age (years)", y = "TERTp Variant Allele\nFrequency (%VAF)",
       color = "Telomere\nlength\nat birth")
p2 <- ggplot(data_plot[mut == "italic(mu) == 1 %*% 10^-6"],
             aes(x = year, y = vaf)) +
  geom_line(aes(group = trial, color = birth_tl), linewidth = 0.3) +
  theme_bw() +
  xlim(0, 90) +
  ylim(0, 0.00413) +
  facet_wrap(~ mut, scales = "free", labeller = label_parsed) +
  scale_color_viridis_c(option = "plasma") +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "gray90", color = "black"),
        strip.text = element_text(size = 10, margin = margin(t = 2, b = 2)),
        text = element_text(family = "Helvetica"),
        legend.position = "none", axis.title.y=element_blank()) +
  labs(x = "Age (years)", y = "TERTp Variant Allele\nFrequency (%VAF)",
       color = "Telomere\nlength\nat birth")
p3 <- ggplot(data_plot[mut == "italic(mu) == 1 %*% 10^-8"],
       aes(x = year, y = vaf)) +
  geom_line(aes(group = trial, color = birth_tl), linewidth = 0.3) +
  theme_bw() +
  xlim(0, 90) +
  ylim(0, 0.000413) +
  facet_wrap(~ mut, scales = "free", labeller = label_parsed) +
  scale_color_viridis_c(option = "plasma") +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "gray90", color = "black"),
        strip.text = element_text(size = 10, margin = margin(t = 2, b = 2)),
        text = element_text(family = "Helvetica"), axis.title.y=element_blank()) +
  labs(x = "Age (years)", y = "TERTp Variant Allele\nFrequency (%VAF)",
       color = "Telomere\nlength\nat birth")

# Combine plots and save
p1 | p2 | p3
ggsave(filename = "popTL_vaf_8500.pdf", width = 9, height = 2.75)