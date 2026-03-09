#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ape)
  library(ggplot2)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Error: Exactly 4 arguments required: <input_dir> <tree.nwk> <species_map_string> <output_dir>")
}

INPUT_DIR  <- args[1]
TREE_FILE  <- args[2]
MAP_STRING <- args[3]
OUTPUT_DIR <- args[4]

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Parsing Species Mapping Dictionary
# Format like "H=Homo_sapiens,M=Mus_musculus"
map_pairs <- strsplit(MAP_STRING, ",")[[1]]
name_map <- list()
for (pair in map_pairs) {
  kv <- strsplit(pair, "=")[[1]]
  name_map[[kv[1]]] <- kv[2]
}
SPECIES_LIST <- names(name_map)

# Load Expression Data & Calculating zTPM
all_ligands <- data.frame()
all_receptors <- data.frame()
all_stats <- data.frame()

for (species in SPECIES_LIST) {
  l_file <- file.path(INPUT_DIR, paste0(species, "_ligand_expression.csv"))
  r_file <- file.path(INPUT_DIR, paste0(species, "_receptor_expression.csv"))
  s_file <- file.path(INPUT_DIR, paste0(species, "_tpm_stats.csv"))
  
  if (!all(file.exists(c(l_file, r_file, s_file)))) {
    warning(paste("Missing input files for species", species, "- Skipping."))
    next
  }
  
  l_data <- read.csv(l_file) %>% mutate(species = species) %>% rename(interactors = 1) 
  r_data <- read.csv(r_file) %>% mutate(species = species) %>% rename(interactors = 1)
  s_data <- read.csv(s_file)
  
  all_ligands <- bind_rows(all_ligands, l_data)
  all_receptors <- bind_rows(all_receptors, r_data)
  all_stats <- bind_rows(all_stats, s_data)
}

# Calculate zTPM = (log2(TPM + 1) - mu) / sigma
all_ligands <- merge(all_ligands, all_stats, by = "species")
all_receptors <- merge(all_receptors, all_stats, by = "species")

all_ligands$ligand_zTPM <- (log2(all_ligands$TPM + 1) - all_ligands$mu) / all_ligands$sigma
all_receptors$receptor_zTPM <- (log2(all_receptors$TPM + 1) - all_receptors$mu) / all_receptors$sigma

full_zTPM <- inner_join(
  all_ligands %>% select(interactors, species, ligand_zTPM),
  all_receptors %>% select(interactors, species, receptor_zTPM),
  by = c("interactors", "species")
) %>% filter(is.finite(ligand_zTPM) & is.finite(receptor_zTPM))

# Map to full scientific names
full_zTPM <- full_zTPM %>%
  mutate(species_fullname = unlist(name_map[species]))

# Load Phylogenetic Tree & Filtering Data
phylo_tree <- read.tree(TREE_FILE)
tree_species_order <- phylo_tree$tip.label

zTPM_clean <- full_zTPM %>% filter(species_fullname %in% tree_species_order)

zTPM_wide <- zTPM_clean %>%
  pivot_wider(
    id_cols = interactors,
    names_from = species_fullname,
    values_from = c(ligand_zTPM, receptor_zTPM),
    names_glue = "{.value}_{species_fullname}"
  ) %>%
  rename_with(~ gsub("ligand_zTPM_", "ligand_exp_", .x)) %>%
  rename_with(~ gsub("receptor_zTPM_", "receptor_exp_", .x))

# Calculate OLS Regression & PIC
stats_results <- list()

for (i in seq_len(nrow(zTPM_wide))) {
  pair_name <- zTPM_wide$interactors[i]
  
  l_vec <- as.numeric(zTPM_wide[i, paste0("ligand_exp_", tree_species_order)])
  r_vec <- as.numeric(zTPM_wide[i, paste0("receptor_exp_", tree_species_order)])
  names(l_vec) <- names(r_vec) <- tree_species_order
  
  valid <- !is.na(l_vec) & !is.na(r_vec)
  if(sum(valid) < 3) next # Need at least 3 valid data points
  
  l_v <- l_vec[valid]; r_v <- r_vec[valid]
  curr_tree <- drop.tip(phylo_tree, names(which(!valid)))
  
  # Ordinary Least Squares (OLS) Regression
  ols_mod <- lm(r_v ~ l_v) 
  ols_s <- summary(ols_mod)
  
  # PIC Regression (forced through origin)
  l_pic <- pic(l_v, curr_tree)
  r_pic <- pic(r_v, curr_tree)
  pic_mod <- lm(r_pic ~ l_pic + 0)
  pic_s <- summary(pic_mod)
  
  # Extract Node Contrasts
  l_contrasts <- as.list(l_pic)
  names(l_contrasts) <- paste0("ligand_pic_node_", seq_along(l_pic))
  r_contrasts <- as.list(r_pic)
  names(r_contrasts) <- paste0("receptor_pic_node_", seq_along(r_pic))
  
  stats_results[[pair_name]] <- cbind(
    data.frame(
      interactors = pair_name,
      intercept = coef(ols_mod)[1],
      slope = coef(ols_mod)[2],
      p_value = if(nrow(ols_s$coefficients) > 1) ols_s$coefficients[2, 4] else NA,
      r_squared = ols_s$r.squared,
      pic_intercept = 0,
      pic_slope = coef(pic_mod)[1],
      pic_p_value = if(nrow(pic_s$coefficients) > 0) pic_s$coefficients[1, 4] else NA,
      pic_r_squared = pic_s$r.squared,
      stringsAsFactors = FALSE
    ),
    as.data.frame(l_contrasts),
    as.data.frame(r_contrasts)
  )
}

final_stats_df <- bind_rows(stats_results)
final_table <- final_stats_df %>% left_join(zTPM_wide, by = "interactors")

output_csv <- file.path(OUTPUT_DIR, "Escalation_Analysis_Table.csv")
fwrite(final_table, output_csv)
cat(sprintf("Results saved to: %s\n", output_csv))

# Plotting
plot_data <- final_table %>%
  filter(!is.na(pic_p_value) & !is.na(slope)) %>%
  mutate(log10_p_value = -log10(pic_p_value)) %>%
  mutate(
    category = "Uncorrelated",
    category = ifelse(slope > 1 & pic_p_value < 0.05, "Symmetric escalation", category),
    category = ifelse(slope < -1 & pic_p_value < 0.05, "Asymmetric escalation", category),
    label = ifelse(category != "Uncorrelated", interactors, NA)
  )

p_value_threshold <- -log10(0.05)

volcano_plot <- ggplot(plot_data, aes(x = slope, y = log10_p_value)) +
  geom_rect(aes(xmin = -Inf, xmax = -1, ymin = -Inf, ymax = Inf), fill = "peachpuff", alpha = 0.05) +
  geom_rect(aes(xmin = 1, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "palegreen", alpha = 0.05) +
  geom_point(aes(color = category), size = 2) +
  scale_color_manual(values = c(
    "Asymmetric escalation" = "red",
    "Symmetric escalation" = "red",
    "Uncorrelated" = "grey50"
  ), name = "Category") +
  geom_hline(yintercept = p_value_threshold, linetype = "dashed", color = "blue") +
  geom_text_repel(aes(label = label), 
                  max.overlaps = Inf, 
                  size = 2.5, 
                  angle = 90, 
                  box.padding = 0.5, 
                  point.padding = 0.2,
                  min.segment.length = 0) +
  labs(
    x = expression("Regression slope (symmetric " * log[10] * ")"),
    y = expression("-log"[10] * "(PIC P value)"),
    title = "Ligand-Receptor Co-evolution Analysis"
  ) +
  annotate("text", x = -4, y = max(plot_data$log10_p_value, na.rm = TRUE) * 0.9, 
           label = "Asymmetric escalation\n(~ correlation)", color = "darkorange", fontface = "bold") +
  annotate("text", x = 2.5, y = max(plot_data$log10_p_value, na.rm = TRUE) * 0.9, 
           label = "Symmetric escalation\n(+ correlation)", color = "darkgreen", fontface = "bold") +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

plot_file <- file.path(OUTPUT_DIR, "coevolution_volcano_plot.pdf")
ggsave(plot_file, plot = volcano_plot, width = 8, height = 6)
cat(sprintf("Volcano plot generated and saved to: %s\n", plot_file))
