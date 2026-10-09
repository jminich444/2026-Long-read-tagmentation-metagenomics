setwd("~/Work/Minich_Collaboration/R_Work/BE_Data/Figure_4_Plots")

library(tidyverse)
library(ggpubr)
library(pals)
library(ggh4x)  # for guide_axis_minor() on 4F

set.seed(42)  # geon-jitter reproducibility

## theme
base_size <- 18
axis_text_size <- 13
viridis_single_color <- viridisLite::viridis(1, begin = 0.5)
text_theme <- theme(
  axis.title   = element_text(size = base_size),
  axis.text    = element_text(size = axis_text_size),
  axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1),
  legend.title = element_text(size = base_size * 0.9),
  legend.text  = element_text(size = base_size * 0.8)
)
plot_theme      <- theme_classic(base_size = base_size) + text_theme  #no gridlines
plot_theme_grid <- theme_minimal(base_size = base_size) + text_theme  # with gridlines
title_theme <- theme(plot.title = element_text(face = "bold", size = base_size * 1.1))

# Saves a titled and an untitled version of the same plot--for use with all plots. 
save_both <- function(plot, filename, title_text, width = 10, height = 6.5, width_no_title = width) {
  ggsave(paste0(filename, "_no_title.png"), plot, width = width_no_title, height = height, dpi = 300)
  ggsave(paste0(filename, "_titled.png"), plot + labs(title = title_text) + title_theme, width = width, height = height, dpi = 300)
}

## Shared theme
final_base_size <- 20   # axis titles and legend titles
final_axis_text <- 18   # axis tick labels
final_legend_text <- 17 # legend text (dropped 1pt from axis_text for legend readability at compact sizes)
final_plot_theme <- theme(
  axis.title       = element_text(size = final_base_size),
  axis.text        = element_text(size = final_axis_text),
  legend.title     = element_text(size = final_base_size),
  legend.text      = element_text(size = final_legend_text),
  legend.key.size  = unit(0.5, "cm")
)

per_sample <- read_csv("~/Work/Minich_Collaboration/Excel Files/Built Env/csv/BE_analysis_sheets_per_sample.csv", show_col_types = FALSE)
per_run    <- read_csv("~/Work/Minich_Collaboration/Excel Files/Built Env/csv/BE_analysis_sheets_per_run.csv", show_col_types = FALSE)

# Base filtered dataset. Drop upstream max_read_length (recomputed from per_run below), add total_gb and mag_category
per_sample_dorm <- per_sample |>
  filter(sample_type == "dorm_sample") |>
  select(-any_of("max_read_length")) |>
  mutate(
    total_gb = number_of_bases / 1e9,
    mag_category = case_when(
      n_complete_MAGs == 0 ~ "0",
      n_complete_MAGs == 1 ~ "1",
      TRUE ~ "Other"
    ),
    mag_category = factor(mag_category, levels = c("0", "1", "Other"))
  )

input_scale <- scale_x_continuous(breaks = seq(0, 280, by = 10))
yield_scale <- scale_y_continuous(limits = c(0, NA), breaks = seq(0, 14, by = 1))


## PLOTS (BE Analysis)
## ____________________
##


dna_input_breaks <- c(0, 15, 30, 45, 60, 75, 150, 300)
dna_input_labels <- c("1-15 ng", "15-30 ng", "30-45 ng", "45-60 ng", "60-75 ng", "75-150 ng", "150-300 ng")

# Add DNA-input bin column so 4E can color points by input bracket
per_sample_dorm_with_bins <- per_sample_dorm |>
  mutate(dna_input_bin = cut(input_DNA_ng, breaks = dna_input_breaks, labels = dna_input_labels, include.lowest = TRUE))

# Per-bin sample counts, used to build "1-15 ng (n=11)" legend labels
bin_counts <- per_sample_dorm_with_bins |> count(dna_input_bin)
dna_input_labels_with_n <- paste0(bin_counts$dna_input_bin, " (n=", bin_counts$n, ")")
names(dna_input_labels_with_n) <- as.character(bin_counts$dna_input_bin)

# Relabel bin factor with per-bin counts embedded
per_sample_dorm_binned <- per_sample_dorm_with_bins |>
  mutate(dna_input_bin = factor(dna_input_labels_with_n[as.character(dna_input_bin)], levels = dna_input_labels_with_n))

# Linear model for predicting yield required to reach a target MAG count.
mag_count_model <- lm(n_complete_MAGs ~ total_gb, data = per_sample_dorm_binned)
print(summary(mag_count_model))

mag_eq <- sprintf(
  "Complete MAGs =\n%.2f + %.2f × yield (Gb)",
  coef(mag_count_model)[1],
  coef(mag_count_model)[2]
)

scale_28 <- list(
  scale_x_continuous(limits = c(0, 28), breaks = seq(0, 28, by = 4)),  # sparser for 6" width
  scale_y_continuous(limits = c(0, 28), breaks = seq(0, 28, by = 4))
)

# 2B - Complete MAGs by sequencing yield
plot_final_4E <- ggplot(per_sample_dorm_binned, aes(x = total_gb, y = n_complete_MAGs)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40", fill = "grey80", linewidth = 0.8) +
  geom_point(aes(color = dna_input_bin), size = 3.5, alpha = 0.85) +
  scale_color_viridis_d(name = "DNA input") +
  stat_cor(method = "spearman", cor.coef.name = "rho",
           label.x.npc = "left", label.y.npc = "top",
           size = final_axis_text * 0.28) +
  annotate("text", x = 27.5, y = 25, label = mag_eq,
           hjust = 1, size = final_axis_text * 0.28) +
  scale_28 +
  labs(x = "Sequencing yield (Gb)", y = "Number of complete MAGs") +
  plot_theme +
  final_plot_theme +
  theme(
    axis.text.x = element_text(angle = 0),  # override 90-deg from base text_theme
    # Legend inside plot area at lower-right 
    legend.position = c(0.98, 0.02),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = alpha("white", 0.85), color = NA)
  )

save_both(plot_final_4E, "final_plots_figure_4E",
          "Complete MAGs by sequencing yield",
          width = 6, height = 6)


# 2C - Sequencing yield by DNA input
# Zero-MAG samples highlighted in red-orange color
# (plotted first in a color defined by the theme, and then zero-MAG points are redrawn on top in a fixed color.)
one_color <- scales::colour_ramp(viridisLite::viridis(256))(1 / max(per_sample_dorm_binned$n_complete_MAGs))
zero_color <- "#C03830"

  input_scale <- scale_x_continuous(breaks = seq(0, 280, by = 40))  # sparser breaks for 6" width, overrides input_scale in the header

plot_final_4D <- ggplot(per_sample_dorm_binned, aes(x = input_DNA_ng, y = total_gb)) +
  geom_point(aes(color = n_complete_MAGs, shape = mag_category), size = 3.5, alpha = 0.85, stroke = 1.2) +
  geom_point(data = filter(per_sample_dorm_binned, mag_category == "0"), aes(shape = mag_category),
             color = zero_color, size = 3.5, alpha = 0.85, stroke = 1.2) +
  scale_color_viridis_c(name = "Complete\nMAGs", breaks = c(0, 5, 10, 15, 20, 25),
                        labels = c("1", "5", "10", "15", "20", "25")) +
  scale_shape_manual(name = NULL, values = c("0" = "\u25B2", "1" = "\u25B2", "Other" = "\u25CF")) +
  guides(
    color = guide_colorbar(order = 1),
    shape = guide_legend(order = 2, override.aes = list(color = c(zero_color, one_color, "black")))
  ) +
  stat_cor(method = "spearman", cor.coef.name = "rho",
           label.x.npc = "left", label.y.npc = "top",
           size = final_axis_text * 0.28) +
  input_scale + yield_scale +
  labs(x = "DNA input (ng)", y = "Sequencing yield (Gb)") +
  plot_theme +
  final_plot_theme +
  theme(
    axis.text.x = element_text(angle = 0),  # override 90-deg from base
    
    legend.position = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = alpha("white", 0.85), color = NA)
  )

save_both(plot_final_4D, "final_plots_figure_4D",
          "Sequencing yield by DNA input",
          width = 6, height = 6)


# Supp. 7B - Read length by DNA input
# Mean and max read lengths overlaid on one plot with per-metric loess fits

# Parse longest_read (integer) from per_run's "longest_read_(with_Q):1" string column
per_run_with_longest <- per_run |>
  mutate(longest_read = as.numeric(str_extract(`longest_read_(with_Q):1`, "^\\d+")))

# Per-sample max read length across all runs of that sample
max_read_length_per_sample <- per_run_with_longest |>
  group_by(sample_id, filename_ID) |>
  summarise(max_read_length = max(longest_read, na.rm = TRUE), .groups = "drop")

# Join per-sample max onto binned dorm data; rename Argo's avg column for brevity
per_sample_dorm_read_lengths <- per_sample_dorm_binned |>
  left_join(max_read_length_per_sample, by = c("sample_id", "filename_ID"))

# Long format for facet_wrap on metric type
per_sample_dorm_read_lengths_long <- per_sample_dorm_read_lengths |>
  select(input_DNA_ng, mean_read_length, max_read_length, n50) |>
  pivot_longer(cols = c(mean_read_length, max_read_length, n50), names_to = "metric", values_to = "read_length") |>
  mutate(metric = factor(case_when(
    metric == "max_read_length"  ~ "Max",
    metric == "n50"              ~ "N50",
    metric == "mean_read_length" ~ "Mean"
  ), levels = c("Max", "N50", "Mean")))

max_cor  <- cor.test(per_sample_dorm_read_lengths$input_DNA_ng, per_sample_dorm_read_lengths$max_read_length,  method = "spearman")
n50_cor  <- cor.test(per_sample_dorm_read_lengths$input_DNA_ng, per_sample_dorm_read_lengths$n50,              method = "spearman")
mean_cor <- cor.test(per_sample_dorm_read_lengths$input_DNA_ng, per_sample_dorm_read_lengths$mean_read_length, method = "spearman")
read_length_colors <- c(
  Max  = viridisLite::viridis(1, begin = 0.15),
  N50  = viridisLite::viridis(1, begin = 0.45),
  Mean = "#800020"
)

plot_final_4F <- ggplot(per_sample_dorm_read_lengths_long, aes(x = input_DNA_ng, y = read_length, color = metric, fill = metric)) +
  geom_hline(yintercept = 5000, linetype = "dotted", color = "grey50", linewidth = 0.4) +
  geom_point(size = 3.5, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.8, alpha = 0.15) +
  scale_color_manual(name = "Read length", values = read_length_colors) +
  scale_fill_manual(name = "Read length", values = read_length_colors) +
  annotate("text", x = 3, y = 58500, hjust = 0, color = read_length_colors[["Max"]],
           size = final_axis_text * 0.28,
           parse = TRUE,
           label = sprintf('rho*" = %.2f, "*italic(p)*" = %.1e"', max_cor$estimate, max_cor$p.value)) +
  annotate("text", x = 3, y = 55500, hjust = 0, color = read_length_colors[["N50"]],
           size = final_axis_text * 0.28,
           parse = TRUE,
           label = sprintf('rho*" = %.2f, "*italic(p)*" = %.1e"', n50_cor$estimate, n50_cor$p.value)) +
  annotate("text", x = 3, y = 52500, hjust = 0, color = read_length_colors[["Mean"]],
           size = final_axis_text * 0.28,
           parse = TRUE,
           label = sprintf('rho*" = %.2f, "*italic(p)*" = %.1e"', mean_cor$estimate, mean_cor$p.value)) +
  scale_x_continuous(breaks = seq(0, 280, by = 40)) +  # sparser breaks for 6" width
  scale_y_continuous(limits = c(0, NA), breaks = seq(0, 60000, by = 5000),
                     minor_breaks = seq(0, 60000, by = 1000),
                     guide = ggh4x::guide_axis_minor()) +
  labs(x = "DNA input (ng)", y = "Read length (bp)") +
  plot_theme +
  final_plot_theme +
  theme(
    axis.text.x = element_text(angle = 0),  
    legend.position = c(0.98, 0.4),
    legend.justification = c("right", "center"),
    legend.background = element_rect(fill = alpha("white", 0.85), color = NA)
  )

save_both(plot_final_4F, "final_plots_figure_4F",
          "Read length by DNA input",
          width = 6, height = 9)


# 2E - Novelty by family

mags <- read_csv("~/Work/Minich_Collaboration/Excel Files/Built Env/csv/BE_analysis_sheets_mags.csv", show_col_types = FALSE)

complete_mags <- mags |>
  filter(Complete == 1, `5S_rRNA` >= 1, `16S_rRNA` >= 1, `23S_rRNA` >= 1, sample_type == "dorm_sample") |>
  mutate(
    family        = str_extract(classification, "f__([^;]*)") |> str_remove("f__"),
    genus         = str_extract(classification, "g__([^;]*)") |> str_remove("g__"),
    species       = str_extract(classification, "s__(.*)$") |> str_remove("s__"),
    novel_species = species == "",
    novel_genus   = genus == ""
  )

by_family_table <- complete_mags |>
  group_by(family) |>
  summarise(
    n_MAGs           = n(),
    n_novel_species   = sum(novel_species),
    pct_novel_species = round(100 * mean(novel_species), 1),
    n_novel_genus     = sum(novel_genus),
    pct_novel_genus   = round(100 * mean(novel_genus), 1),
    .groups = "drop"
  ) |>
  arrange(desc(n_MAGs))

family_label_size <- 3.2

plot_final_4G <- ggplot(by_family_table, aes(x = n_MAGs, y = pct_novel_species)) +
  geom_point(aes(size = n_novel_genus, color = n_novel_genus), alpha = 0.8) +
  scale_color_viridis_c(name = "Genus-novel\nMAGs", breaks = sort(unique(by_family_table$n_novel_genus))) +
  scale_size_continuous(name = "Genus-novel\nMAGs", range = c(2, 10), breaks = sort(unique(by_family_table$n_novel_genus))) +
  guides(color = guide_legend(), size = guide_legend()) +
  ggrepel::geom_text_repel(
    aes(label = family), size = family_label_size, max.overlaps = Inf,
    force = 10, force_pull = 0.5, box.padding = 1.0, point.padding = 0.3,
    segment.size = 0.25, segment.color = "grey60", segment.alpha = 0.7,
    min.segment.length = 0, seed = 1
  ) +
  scale_y_continuous(limits = c(-14.4, 114.4), breaks = seq(0, 100, by = 25)) +
  labs(x = "Number of complete MAGs in family", y = "Novel at species level (%)") +
  plot_theme +
  final_plot_theme

save_both(plot_final_4G, "final_plots_figure_4G",
          "Novelty by family",
          width = 11, height = 7.29)




# SUPPLEMENTARY PLOTS
# ________________________

# s7A - Percentage of human DNA removed
plot_final_4_s6D <- ggplot(per_run |> filter(sample_type == "dorm_sample"), aes(x = "", y = `Percentage removed`)) +
  geom_boxplot(width = 0.06, fill = "grey70", outlier.shape = NA) +
  geom_jitter(width = 0.04, size = 1.8, alpha = 0.7) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.01), breaks = c(0, 0.01, 0.1, 1, 6)) +
  labs(x = NULL, y = "Human DNA removed (%, pseudo-log scale)") +
  plot_theme +
  theme(axis.text.x = element_blank())

save_both(plot_final_4_s6D, "final_plots_figure_4_s6D", "Percentage of human DNA removed", width = 7, height = 7, width_no_title = 4)


# s7E - Novelty by genus (bubble)
by_genus_table <- complete_mags |>
  filter(!novel_genus) |>
  group_by(genus) |>
  summarise(
    n_MAGs            = n(),
    n_novel_species    = sum(novel_species),
    pct_novel_species  = round(100 * mean(novel_species), 1),
    .groups = "drop"
  )

by_genus_coords <- by_genus_table |>
  group_by(n_MAGs, pct_novel_species) |>
  summarise(
    n_genera_here = n(),
    label         = str_wrap(paste(sort(genus), collapse = ", "), width = 45),
    .groups = "drop"
  )

plot_final_4_s6E <- ggplot(by_genus_coords, aes(x = n_MAGs, y = pct_novel_species)) +
  geom_point(aes(size = n_genera_here), color = viridisLite::viridis(1, begin = 0.5), alpha = 0.8) +
  scale_size_continuous(name = "Genera at\nthis point", range = c(2, 9), breaks = c(1, 5, 10, 15, 20, 26)) +
  ggrepel::geom_text_repel(
    aes(label = label), size = family_label_size, max.overlaps = Inf, lineheight = 0.9,
    force = 60, force_pull = 0.3, box.padding = 1.2, point.padding = 0.3,
    segment.size = 0.25, segment.color = "grey60", segment.alpha = 0.7,
    min.segment.length = 0, seed = 1
  ) +
  labs(x = "Number of complete MAGs in genus", y = "Novel at species level (%)") +
  plot_theme

save_both(plot_final_4_s6E, "final_plots_figure_4_s6E", "Novelty by genus", width = 14, height = 10)


# s7C - GC content
plot_final_4_s6A <- ggplot(complete_mags, aes(x = "", y = GC_Content * 100)) +
  geom_boxplot(width = 0.06, fill = "grey70", outlier.shape = NA) +
  geom_jitter(width = 0.04, size = 1.8, alpha = 0.7) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  labs(x = NULL, y = "GC content (%)") +
  plot_theme +
  theme(axis.text.x = element_blank())

save_both(plot_final_4_s6A, "final_plots_figure_4_s6A", "GC content", width = 7, height = 7, width_no_title = 4)

# s7D - Genome size
plot_final_4_s6B <- ggplot(complete_mags, aes(x = "", y = Genome_Size / 1e6)) +
  geom_boxplot(width = 0.06, fill = "grey70", outlier.shape = NA) +
  geom_jitter(width = 0.04, size = 1.8, alpha = 0.7) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  labs(x = NULL, y = "Genome size (Mb)") +
  plot_theme +
  theme(axis.text.x = element_blank())

save_both(plot_final_4_s6B, "final_plots_figure_4_s6B", "Genome size", width = 7, height = 7, width_no_title = 4)
