setwd("~/Work/Minich_Collaboration/R_Work/BE_Data/Figure_4_Plots")

library(tidyverse)
library(ggpubr)
library(pals)
library(ggh4x)  # for force_panelsizes() to make plots share aspect ratio

# Shared theme: large text for two-panel layouts, dense axis ticks. 
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
plot_theme      <- theme_classic(base_size = base_size) + text_theme  # no gridlines
plot_theme_grid <- theme_minimal(base_size = base_size) + text_theme  # with gridlines
title_theme <- theme(plot.title = element_text(face = "bold", size = base_size * 1.1))

# Saves a titled and an untitled version of the same plot. 

save_both <- function(plot, filename, title_text, width = 10, height = 6.5, width_no_title = width) {
  ggsave(paste0(filename, "_no_title.png"), plot, width = width_no_title, height = height, dpi = 300)
  ggsave(paste0(filename, "_titled.png"), plot + labs(title = title_text) + title_theme, width = width, height = height, dpi = 300)
}


# Data
argo_community_profile <- read_csv("~/Work/Minich_Collaboration/Excel Files/Built Env/csv/BE_analysis_sheets_Argo_community_profile.csv", show_col_types = FALSE)
argo_sarg_amr_genes    <- read_csv("~/Work/Minich_Collaboration/Excel Files/Built Env/csv/BE_analysis_sheets_Argo_sarg_amr_genes.csv", show_col_types = FALSE)
mags_fp                <- read_csv("~/Work/Minich_Collaboration/Excel Files/Built Env/csv/BE_analysis_sheets_mags.csv", show_col_types = FALSE) |>
  filter(sample_type != "control")  # exclude control MAGs
valid_mag_ids_fp       <- paste0(mags_fp$Sample, ".", mags_fp$Name)
reads_mag_arg_all_fp   <- read_csv("~/Work/Minich_Collaboration/Excel Files/Built Env/csv/BE_analysis_sheets_Argo_read_mapping.csv", show_col_types = FALSE) |>
  filter(!filename_ID %in% c("H07", "G07"),  # exclude reads from H07/G07 controls
         mag_id %in% valid_mag_ids_fp)  # exclude reads mapped to control MAGs
reads_mag_arg_fp       <- reads_mag_arg_all_fp |> filter(read_classification_status == "unclassified")



## Shared setup 
argo_sarg_dorm_fp <- argo_sarg_amr_genes |>
  filter(sample_type == "dorm_sample") |>
  mutate(species_clean = if_else(lineage == "unclassified", "Unclassified", str_extract(lineage, "[^;]+$")))

## ARG type ordering (ascending total copy abundance)
arg_type_order_fp <- argo_sarg_dorm_fp |>
  group_by(type) |>
  summarise(total_copy = sum(copy), .groups = "drop") |>
  arrange(total_copy) |>
  pull(type)


## Shared setup for s6F 
type_species_rel_copy_fp <- argo_sarg_dorm_fp |>
  group_by(type, species_clean) |>
  summarise(copy_sum = sum(copy), .groups = "drop") |>
  group_by(type) |>
  mutate(rel_copy = copy_sum / sum(copy_sum)) |>
  ungroup()

top_species_copy_fp <- argo_sarg_dorm_fp |>
  filter(species_clean != "Unclassified") |>
  group_by(species_clean) |>
  summarise(total_copy = sum(copy), .groups = "drop") |>
  slice_max(total_copy, n = 20) |>
  pull(species_clean)

species_bucket_levels_copy_fp <- c(top_species_copy_fp, "Other", "Unclassified")

type_species_composition_copy_fp <- type_species_rel_copy_fp |>
  mutate(species_bucket = case_when(
    species_clean == "Unclassified" ~ "Unclassified",
    species_clean %in% top_species_copy_fp ~ species_clean,
    TRUE ~ "Other"
  )) |>
  group_by(type, species_bucket) |>
  summarise(rel_copy = sum(rel_copy), .groups = "drop") |>
  mutate(
    species_bucket = factor(species_bucket, levels = species_bucket_levels_copy_fp),
    type = factor(type, levels = arg_type_order_fp)
  )

species_colors_copy_fp <- c(as.vector(stepped2(length(top_species_copy_fp))), "grey85", "grey60")
names(species_colors_copy_fp) <- species_bucket_levels_copy_fp

## shared font sizes
final_base_size <- 20   ## axis titles and legend titles
final_axis_text <- 18   # axis tick labels and legend text

final_plot_theme <- theme(
  axis.title       = element_text(size = final_base_size),
  axis.text        = element_text(size = final_axis_text),
  legend.title     = element_text(size = final_base_size),
  legend.text      = element_text(size = final_axis_text),
  legend.key.size  = unit(0.5, "cm")
)


# mean cpg and top-20 species by summed load. 
sarg_species_argtype_mean_cpg <- argo_sarg_dorm_fp |>
  # Detected only: rows in Argo_sarg_amr_genes are already all detected
  # events (the sheet doesn't contain zero-cpg rows), so this is a
  #just a safety filter.
  filter(abundance > 0) |>
  group_by(species_clean, subtype, type) |>
  # Step 2: mean cpg across samples where this species has this subtype
  summarise(mean_cpg_subtype = mean(abundance), .groups = "drop") |>
  group_by(species_clean, type) |>
  # Step 3: sum subtype mean-cpg within each (species, drug class)
  summarise(species_type_load = sum(mean_cpg_subtype), .groups = "drop")

# Rank species by total mean-cpg across all classes (top 20 by
# per-genome resistance load, per-genome-normalized weighting)
top_species_5_06 <- sarg_species_argtype_mean_cpg |>
  filter(species_clean != "Unclassified") |>
  group_by(species_clean) |>
  summarise(total_load = sum(species_type_load), .groups = "drop") |>
  slice_max(total_load, n = 20) |>
  arrange(desc(total_load)) |>
  pull(species_clean)



###### ------- Plots---------

# Supp. 7G - Species-level community composition per sample

community_dorm_5_07 <- argo_community_profile |>
  filter(sample_type == "dorm_sample") |>
  mutate(
    species_clean = if_else(
      str_detect(species, "nclassified") | is.na(species) | species == "",
      "Unclassified",
      species
    ),
    Dorm = factor(Dorm, levels = c("Dawson", "Memorial", "Collins", "NoRo",
                                   "Penland", "SoRo", "Earl", "Teal", "Heritage"))
  )

top_species_5_07 <- community_dorm_5_07 |>
  filter(species_clean != "Unclassified") |>
  group_by(species_clean) |>
  summarise(mean_abundance = mean(abundance), .groups = "drop") |>
  slice_max(mean_abundance, n = 20) |>
  pull(species_clean)

species_bucket_levels_5_07 <- c(top_species_5_07, "Other", "Unclassified")

species_composition_5_07 <- community_dorm_5_07 |>
  mutate(species_bucket = case_when(
    species_clean == "Unclassified" ~ "Unclassified",
    species_clean %in% top_species_5_07 ~ species_clean,
    TRUE ~ "Other"
  )) |>
  group_by(sample_id, Dorm, species_bucket) |>
  summarise(abundance = sum(abundance), .groups = "drop") |>
  mutate(species_bucket = factor(species_bucket, levels = species_bucket_levels_5_07))

sample_order_5_07 <- species_composition_5_07 |>
  distinct(sample_id, Dorm) |>
  arrange(Dorm, sample_id) |>
  pull(sample_id)

species_composition_5_07 <- species_composition_5_07 |>
  mutate(sample_id = factor(sample_id, levels = sample_order_5_07))

species_colors_5_07 <- c(as.vector(stepped2(length(top_species_5_07))), "grey85", "grey60")
names(species_colors_5_07) <- species_bucket_levels_5_07

plot_final_5A <- ggplot(species_composition_5_07,
                          aes(x = sample_id, y = abundance, fill = species_bucket)) +
  geom_col(width = 0.85) +
  facet_grid(. ~ Dorm, scales = "free_x", space = "free_x") +
  scale_fill_manual(name = "Species", values = species_colors_5_07) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Samples (grouped by Dorm, oldest to youngest)", y = "Relative abundance") +
  plot_theme +
  final_plot_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.text = element_text(size = final_axis_text, face = "italic"),
    strip.text = element_text(size = final_axis_text * 0.75, face = "bold", angle = 90),
    panel.spacing = unit(0.3, "lines")
  )

save_both(plot_final_5A, "final_plots_figure_5A",
          "Species-level community composition per sample",
          width = 14, height = 8)


# Supp. 7J - Relative per-genome ARG load by drug class

# Species presence from Argo community profile
species_presence_5_10 <- argo_community_profile |>
  filter(sample_type == "dorm_sample", abundance > 0) |>
  mutate(species_clean = if_else(
    str_detect(species, "nclassified") | is.na(species) | species == "",
    "Unclassified", species
  )) |>
  filter(species_clean != "Unclassified") |>
  distinct(species_clean, sample_id)

# Per-sample per-class cpg burden for detected (species, subtype, sample) events
per_sample_class_cpg_5_10 <- argo_sarg_dorm_fp |>
  filter(species_clean != "Unclassified") |>
  group_by(species_clean, sample_id, type) |>
  summarise(class_cpg_sum = sum(abundance), .groups = "drop")

# Complete grid; for every (species, presence-sample) join every ARG
# class, filling missing (present but no detection) combinations with 0; this means we include all samples where the organism was found, even if no resistance in every sample
all_types_5_10 <- unique(argo_sarg_dorm_fp$type)

species_class_5_10 <- species_presence_5_10 |>
  crossing(type = all_types_5_10) |>
  left_join(per_sample_class_cpg_5_10,
            by = c("species_clean", "sample_id", "type")) |>
  mutate(class_cpg_sum = replace_na(class_cpg_sum, 0))

# Mean per (species, class) across presence samples
species_class_burden_5_10 <- species_class_5_10 |>
  group_by(species_clean, type) |>
  summarise(mean_burden = mean(class_cpg_sum), .groups = "drop")

# Overall species ranking (sum mean_burden across classes)
top_species_5_10 <- species_class_burden_5_10 |>
  group_by(species_clean) |>
  summarise(total_burden = sum(mean_burden), .groups = "drop") |>
  slice_max(total_burden, n = 20) |>
  arrange(desc(total_burden)) |>
  pull(species_clean)

species_bucket_levels_5_10 <- c(top_species_5_10, "Other")

# Pool non-top species into "Other" and compute class denominators
class_species_composition_5_10 <- species_class_burden_5_10 |>
  mutate(species_bucket = if_else(
    species_clean %in% top_species_5_10, species_clean, "Other"
  )) |>
  group_by(type, species_bucket) |>
  summarise(bucket_burden = sum(mean_burden), .groups = "drop") |>
  group_by(type) |>
  mutate(rel_burden = bucket_burden / sum(bucket_burden)) |>
  ungroup() |>
  mutate(
    species_bucket = factor(species_bucket, levels = species_bucket_levels_5_10),
    type = factor(type, levels = arg_type_order_fp)
  )

species_colors_5_10 <- c(as.vector(stepped2(length(top_species_5_10))), "grey85")
names(species_colors_5_10) <- species_bucket_levels_5_10

plot_final_5B <- ggplot(class_species_composition_5_10,
                          aes(x = rel_burden, y = type, fill = species_bucket)) +
  geom_col(width = 0.7) +
  scale_fill_manual(name = "Species", values = species_colors_5_10) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Relative per-genome ARG load", y = NULL) +
  plot_theme +
  final_plot_theme +
  theme(
    axis.text.x = element_text(angle = 0),
    legend.text = element_text(size = final_axis_text, face = "italic")
  )

save_both(plot_final_5B, "final_plots_figure_5B",
          "Relative per-genome ARG load",
          width = 14, height = 8)





# Supp. 7K - Chromosomal and plasmid ARGs by Argo classification, with paired Wilcoxon tests
carrier_by_class_fp <- argo_sarg_dorm_fp |>
  mutate(classification_status = if_else(lineage == "unclassified", "Unclassified", "Classified")) |>
  group_by(classification_status, carrier) |>
  summarise(n = round(sum(copy)), .groups = "drop") |>
  group_by(classification_status) |>
  mutate(prop = n / sum(n)) |>
  ungroup()

# Statistical section: paired Wilcoxon signed-rank on per-sample copy
# sums, one comparison per carrier type. See if differences in chromosome and plasmid counts are significant

per_sample_5_04 <- argo_sarg_dorm_fp |>
  mutate(classification_status = if_else(lineage == "unclassified", "Unclassified", "Classified")) |>
  group_by(sample_id, classification_status, carrier) |>
  summarise(copies = sum(copy), .groups = "drop") |>
  complete(sample_id, classification_status, carrier, fill = list(copies = 0)) |>
  mutate(group = paste(classification_status, carrier, sep = "_")) |>
  select(sample_id, group, copies) |>
  pivot_wider(names_from = group, values_from = copies, values_fill = 0)

wilcox_chr_5_04 <- wilcox.test(per_sample_5_04$Classified_chromosome,
                                per_sample_5_04$Unclassified_chromosome,
                                paired = TRUE, exact = FALSE)
wilcox_pla_5_04 <- wilcox.test(per_sample_5_04$Classified_plasmid,
                                per_sample_5_04$Unclassified_plasmid,
                                paired = TRUE, exact = FALSE)

# Significance-star notation 
sig_stars <- function(p) {
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
}
stars_chr_5_04 <- sig_stars(wilcox_chr_5_04$p.value)
stars_pla_5_04 <- sig_stars(wilcox_pla_5_04$p.value)

# Code for brackets to be displayed in non-typical locations. Two rows above the top of the plot area. Chromosome
# comparison is drawn higher (y = 1.18); plasmid comparison sits just
# below it (y = 1.08). 

bracket_5_04 <- tibble(
  y_bracket = c(1.08, 1.18),
  y_label   = c(1.11, 1.21),
  label     = c(paste0("plasmid ", stars_pla_5_04),
                paste0("chromosome ", stars_chr_5_04)),
  color     = c("#DD8452", "#4C72B0")
)

plot_final_5E <- ggplot(carrier_by_class_fp, aes(x = classification_status, y = prop, fill = carrier)) +
  geom_col(width = 0.87) +
  geom_text(data = ~ filter(.x, carrier == "chromosome"),
            aes(label = paste0(carrier, "\n", scales::percent(prop, accuracy = 0.1), "\nn=", n)),
            position = position_stack(vjust = 0.5), color = "white",
            size = final_axis_text * 0.25) +
  geom_text(data = ~ filter(.x, carrier == "plasmid"),
            aes(y = 0.09,  # fixed anchor above the 0-baseline; higher than mid-segment position so all three lines of the label (plasmid / X% / n=Y) fit above the x-axis without clipping
                label = paste0(carrier, "\n", scales::percent(prop, accuracy = 0.1), "\nn=", n)),
            color = "white", size = final_axis_text * 0.25) +
  # Significance brackets. 
  geom_segment(data = bracket_5_04, inherit.aes = FALSE,
               aes(x = 1, xend = 2, y = y_bracket, yend = y_bracket, color = color),
               linewidth = 0.6, show.legend = FALSE) +
  geom_segment(data = bracket_5_04, inherit.aes = FALSE,
               aes(x = 1, xend = 1, y = y_bracket - 0.02, yend = y_bracket, color = color),
               linewidth = 0.6, show.legend = FALSE) +
  geom_segment(data = bracket_5_04, inherit.aes = FALSE,
               aes(x = 2, xend = 2, y = y_bracket - 0.02, yend = y_bracket, color = color),
               linewidth = 0.6, show.legend = FALSE) +
  geom_text(data = bracket_5_04, inherit.aes = FALSE,
            aes(x = 1.5, y = y_label, label = label, color = color),
            size = final_axis_text * 0.28, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_manual(values = c(chromosome = "#4C72B0", plasmid = "#DD8452"), guide = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                     expand = expansion(mult = c(0, 0.28)),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     oob = scales::oob_keep) +
  labs(x = NULL, y = "Proportion of ARG copies") +
  plot_theme +
  final_plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  
  # Force absolute panel dimensions so 5E and 5F match l
  # when exported at the same figure size
  
  ggh4x::force_panelsizes(rows = unit(6, "in"), cols = unit(2.64, "in"))

save_both(plot_final_5E, "final_plots_figure_5E",
          "Chromosomal and plasmid ARGs by Argo classification status",
          width = 14/3, height = 8)


# 2F - Chromosomal and plasmid assignment of ARG-containing reads by drug class
carrier_by_type_fp <- argo_sarg_dorm_fp |>
  group_by(type, carrier) |>
  summarise(copy_sum = sum(copy), .groups = "drop") |>
  group_by(type) |>
  mutate(
    prop = copy_sum / sum(copy_sum),
    total_copies = sum(copy_sum)
  ) |>
  ungroup() |>
  mutate(type = factor(type, levels = arg_type_order_fp))

# Wide-format labels: one row per type with rounded chromosome and
#plasmid counts as separate columns so that each number can be drawn in its matching bar color 
n_labels_5_03 <- carrier_by_type_fp |>
  mutate(copy_rounded = round(copy_sum)) |>
  select(type, carrier, copy_rounded) |>
  pivot_wider(names_from = carrier, values_from = copy_rounded, values_fill = 0)

plot_final_5D <- ggplot(carrier_by_type_fp, aes(x = type, y = prop, fill = carrier)) +
  geom_col(width = 0.4) +

  geom_text(data = n_labels_5_03,
            aes(x = type, y = 1.13, label = plasmid),
            inherit.aes = FALSE, hjust = 1, color = "#DD8452",
            size = final_axis_text * 0.28) +
  geom_text(data = n_labels_5_03,
            aes(x = type, y = 1.16, label = "/"),
            inherit.aes = FALSE, hjust = 0.5, color = "grey40",
            size = final_axis_text * 0.28) +
  geom_text(data = n_labels_5_03,
            aes(x = type, y = 1.30, label = chromosome),
            inherit.aes = FALSE, hjust = 1, color = "#4C72B0",
            size = final_axis_text * 0.28) +
  coord_flip(clip = "off") +
  scale_fill_manual(name = "Carrier", values = c(chromosome = "#4C72B0", plasmid = "#DD8452")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                     expand = expansion(mult = c(0, 0.36)),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     oob = scales::oob_keep) +
  labs(x = NULL, y = "ARG sources") +
  plot_theme +
  final_plot_theme +
  theme(
    axis.text.x = element_text(angle = 0),
    axis.title.x = element_text(hjust = 0.42)  # visual re-center under the plot area (label area on the right pushes default center rightward)
  )

save_both(plot_final_5D, "final_plots_figure_5D",
          "Chromosomal and plasmid assignment of ARG-containing reads by drug class",
          width = 14, height = 8)


# 5G - ARG-containing reads by post-MAG-mapping classification status
status_order_fp <- c("GTDB-classified MAG", "Species-novel MAG")

reads_mag_arg_fp <- reads_mag_arg_fp |>
  mutate(mag_status = if_else(species_novel, "Species-novel MAG", "GTDB-classified MAG"))

status_summary_fp <- reads_mag_arg_fp |>
  group_by(mag_status) |>
  summarise(n_unique_reads = n_distinct(paste(filename_ID, read_id)), .groups = "drop") |>
  mutate(mag_status = factor(mag_status, levels = status_order_fp))

status_colors_fp <- as.vector(stepped2(20))[c(1, 5)]
names(status_colors_fp) <- status_order_fp

plot_final_5F <- ggplot(status_summary_fp, aes(x = mag_status, y = n_unique_reads, fill = mag_status)) +
  geom_col(width = 0.435 / 0.565) +
  geom_text(aes(label = n_unique_reads), vjust = -0.5, size = final_axis_text * 0.35) +
  scale_fill_manual(values = status_colors_fp, guide = "none") +
  scale_x_discrete(expand = expansion(add = 0.3825 / 0.565)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Unique ARG-containing reads (n)") +
  plot_theme +
  final_plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  # Again, just matches 5E and 5F

  ggh4x::force_panelsizes(rows = unit(6, "in"), cols = unit(2.64 * 1.33 / 2.2, "in"))

save_both(plot_final_5F, "final_plots_figure_5F",
          "ARG-containing reads by post-MAG-mapping classification",
          width = 14/3 - 2.64 * 0.87 / 2.2, height = 8)


# 5Ha - ARG types mapped to GTDB-classified species (heatmap)

organism_lookup_fp <- mags_fp |>
  filter(Complete == 1, `5S_rRNA` > 0, `16S_rRNA` > 0, `23S_rRNA` > 0) |>
  filter(!str_detect(classification, "s__$")) |>
  mutate(
    mag_id = paste0(Sample, ".", Name),
    organism_id = str_remove(str_extract(classification, "s__.*$"), "^s__")
  ) |>
  select(mag_id, organism_id)

organism_by_argtype_all_fp <- reads_mag_arg_all_fp |>
  inner_join(organism_lookup_fp, by = "mag_id") |>
  group_by(organism_id, arg_type) |>
  summarise(
    total_reads = n_distinct(paste(filename_ID, read_id)),
    unclassified_reads = n_distinct(paste(filename_ID, read_id)[read_classification_status == "unclassified"]),
    .groups = "drop"
  ) |>
  mutate(cell_label = sprintf("%d (%d)", total_reads, unclassified_reads))

organism_order_5_09 <- organism_by_argtype_all_fp |>
  group_by(organism_id) |>
  summarise(total_unclassified = sum(unclassified_reads),
            total_all = sum(total_reads),
            .groups = "drop") |>
  arrange(total_unclassified, total_all) |>
  pull(organism_id)

argtype_order_5_09 <- organism_by_argtype_all_fp |>
  group_by(arg_type) |>
  summarise(total_unclassified = sum(unclassified_reads),
            total_all = sum(total_reads),
            .groups = "drop") |>
  arrange(total_unclassified, total_all) |>
  pull(arg_type)

organism_by_argtype_all_fp <- organism_by_argtype_all_fp |>
  mutate(
    organism_id = factor(organism_id, levels = organism_order_5_09),
    arg_type = factor(arg_type, levels = argtype_order_5_09)
  )

plot_final_5G <- ggplot(organism_by_argtype_all_fp,
                        aes(x = arg_type, y = organism_id,
                            fill = ifelse(unclassified_reads == 0, NA_real_, unclassified_reads))) +
  geom_tile(color = "white") +
  geom_text(aes(label = cell_label), size = final_axis_text * 0.136) +
  scale_fill_viridis_c(name = "Unclassified\nreads (log10)",
                       trans = "log10", na.value = "grey90") +
  labs(x = "ARG type", y = NULL) +
  plot_theme +
  final_plot_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                               size = final_axis_text - 3),
    axis.text.y = element_text(size = final_axis_text - 3, face = "italic"),
    legend.title = element_text(size = final_axis_text - 2),
    legend.text = element_text(size = final_axis_text - 3),
    legend.position = c(0.02, 0.02),
    legend.justification = c("left", "bottom"),
    legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
    legend.key.size = unit(0.6, "cm")
  )

save_both(plot_final_5G, "final_plots_figure_5G",
          "ARG types mapped to GTDB-classified species (all reads, unclassified recovery highlighted)",
          width = 12, height = 12)




# 5H - ARG types mapped to GTDB-classified species (heatmap)
# Bi-clustered version (Euclidean, average linkage on log10 reads)

# Matrix form: rows = species, columns = ARG types, values = total reads (0 where absent)
mat_5_07b <- organism_by_argtype_all_fp |>
  select(organism_id, arg_type, total_reads) |>
  pivot_wider(names_from = arg_type, values_from = total_reads, values_fill = 0) |>
  column_to_rownames("organism_id") |>
  as.matrix()

# Cluster on log10(x+1) to compress the 4-order-of-magnitude range
mat_5_07b_log <- log10(mat_5_07b + 1)
row_order_5_07b <- rownames(mat_5_07b)[hclust(dist(mat_5_07b_log), method = "average")$order]
col_order_5_07b <- colnames(mat_5_07b)[hclust(dist(t(mat_5_07b_log)), method = "average")$order]

organism_by_argtype_clustered_fp <- organism_by_argtype_all_fp |>
  mutate(
    organism_id = factor(as.character(organism_id), levels = row_order_5_07b),
    arg_type    = factor(as.character(arg_type),    levels = col_order_5_07b)
  )

plot_final_5G_b <- ggplot(organism_by_argtype_clustered_fp,
                           aes(x = arg_type, y = organism_id,
                               fill = ifelse(unclassified_reads == 0, NA_real_, unclassified_reads))) +
  geom_tile(color = "white") +
  geom_text(aes(label = cell_label), size = final_axis_text * 0.136) +
  scale_fill_viridis_c(name = "Unclassified\nreads (log10)",
                       trans = "log10", na.value = "grey90") +
  labs(x = "ARG type", y = NULL) +
  plot_theme +
  final_plot_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                               size = final_axis_text - 3),
    axis.text.y = element_text(size = final_axis_text - 3, face = "italic"),
    legend.position = "none"
  )

save_both(plot_final_5G_b, "final_plots_figure_5G-b",
          "ARG types mapped to GTDB-classified species (bi-clustered)",
          width = 12, height = 12)


# SUPPLEMENTARY PLOTS
# ________________________

# s7F - Total ARG copies by drug class, summed across all species and samples
class_copies_5_11 <- argo_sarg_dorm_fp |>
  group_by(type) |>
  summarise(total_copies = round(sum(copy)), .groups = "drop") |>
  arrange(desc(total_copies)) |>
  mutate(type = factor(type, levels = rev(type)))  # rev so descending shows top-down when horizontal

plot_final_5_s6F <- ggplot(class_copies_5_11, aes(x = total_copies, y = type)) +
  geom_col(width = 0.5, fill = "grey60") +
  geom_text(aes(label = format(total_copies, big.mark = ",")),
            hjust = -0.15, size = final_axis_text * 0.28) +
  scale_x_log10(
    limits = c(1, 3e5),
    breaks = NULL,
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +  # allow labels to overflow the axis line without being clipped
  labs(x = NULL, y = NULL, caption = "Total ARG copies across all species (log10)") +
  plot_theme +
  final_plot_theme +
  theme(
    axis.text.x = element_text(angle = 0),
    plot.caption = element_text(hjust = 0.5, size = final_base_size, margin = margin(t = 10)),
    plot.caption.position = "plot",
    plot.margin = margin(t = 10, r = 60, b = 10, l = 10)
  )

save_both(plot_final_5_s6F, "final_plots_figure_5_s6F",
          "Total ARG copies by drug class across all species",
          width = 14/2, height = 8)


# s7G - Species distribution across ARG types (by copy count)
plot_final_5_s6G <- ggplot(type_species_composition_copy_fp,
                          aes(x = rel_copy, y = type, fill = species_bucket)) +
  geom_col(width = 0.7) +
  scale_fill_manual(name = "Species", values = species_colors_copy_fp) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = "Relative abundance based on copy count", y = NULL) +
  plot_theme +
  final_plot_theme +
  theme(
    axis.text.x = element_text(angle = 0),
    legend.text = element_text(size = final_axis_text, face = "italic")  # italicize species names in legend
  )

save_both(plot_final_5_s6G, "final_plots_figure_5_s6G",
          "Species distribution across ARG types",
          width = 14, height = 8)



# s7I - Percentage of ARG-carrying species resistant to each drug class
arg_carriers_by_class_5_s1 <- argo_sarg_dorm_fp |>
  filter(species_clean != "Unclassified") |>
  group_by(type) |>
  summarise(n_species = n_distinct(species_clean), .groups = "drop") |>
  mutate(pct = 100 * n_species / n_distinct(argo_sarg_dorm_fp$species_clean[argo_sarg_dorm_fp$species_clean != "Unclassified"])) |>
  arrange(desc(pct)) |>
  mutate(type = factor(type, levels = rev(type)))

plot_final_5_s6H <- ggplot(arg_carriers_by_class_5_s1, aes(x = pct, y = type)) +
  geom_col(width = 0.3, fill = "grey60") +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            hjust = -0.15, size = final_axis_text * 0.28) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 50, 100),
                     labels = c("0", "50", "100"), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = "% of classified ARG carriers", y = NULL) +
  plot_theme +
  final_plot_theme +
  theme(
    axis.text.x = element_text(angle = 0),
    plot.margin = margin(t = 10, r = 60, b = 10, l = 10)
  )

save_both(plot_final_5_s6H, "final_plots_figure_5_s6H",
          "Percentage of classified ARG-carrying species resistant to each drug class",
          width = 14/2, height = 8)
