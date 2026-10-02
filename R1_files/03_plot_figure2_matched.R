# 03_plot_figure2_matched.R
# Generates publication-ready Figure 2 heatmap and saves to PDF and PNG.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)
  library(gtable)
})

# 1. Load data from current directory
critria_file <- list.files(pattern = "Supplemental_File_4_criTRia-curations.*\\.tsv")[1]
if (is.na(critria_file)) stop("Supplemental_File_4 criTRia file not found in current directory!")
critria <- read_tsv(critria_file, show_col_types = FALSE)
matrix_wide <- readRDS("discordance_matrix_wide.rds")

# Column order
target_submitters <- c("criTRia", "Labcorp", "Orphanet", "Ambry", "PanelApp", "G2P", "ClinGen", "Lab MM", "Illumina", "Myriad")

# 2. Count evaluated loci per submitter
counts <- sapply(target_submitters, function(grp) {
  sum(matrix_wide[[grp]] != "Uncurated")
})
counts["Conflict"] <- sum(matrix_wide$Conflict)

cat("Evaluated counts per column:\n")
print(counts)

# Create column labels with count annotations on top
col_labels <- c(
  "criTRia" = paste0("criTRia\n(n=", counts["criTRia"], ")"),
  "Labcorp" = paste0("Labcorp\n(n=", counts["Labcorp"], ")"),
  "Orphanet" = paste0("Orphanet\n(n=", counts["Orphanet"], ")"),
  "Ambry" = paste0("Ambry\n(n=", counts["Ambry"], ")"),
  "PanelApp" = paste0("PanelApp\n(n=", counts["PanelApp"], ")"),
  "G2P" = paste0("G2P\n(n=", counts["G2P"], ")"),
  "ClinGen" = paste0("ClinGen\n(n=", counts["ClinGen"], ")"),
  "Lab MM" = paste0("Lab MM\n(n=", counts["Lab MM"], ")"),
  "Illumina" = paste0("Illumina\n(n=", counts["Illumina"], ")"),
  "Myriad" = paste0("Myriad\n(n=", counts["Myriad"], ")"),
  "Conflict" = paste0("Conflict\n(n=", counts["Conflict"], ")")
)

# 3. Format heatmap long dataset
# Ensure loci are sorted alphabetically from A (top) to Z (bottom)
locus_levels_rev <- sort(unique(matrix_wide$Locus_ID), decreasing = TRUE)

heatmap_df <- matrix_wide %>%
  mutate(Conflict_Status = ifelse(Conflict, "Conflict", "No Conflict")) %>%
  select(Locus_ID, all_of(target_submitters), Conflict = Conflict_Status) %>%
  pivot_longer(cols = c(all_of(target_submitters), "Conflict"), names_to = "Group", values_to = "Score") %>%
  mutate(
    Group_Label = factor(col_labels[Group], levels = unname(col_labels)),
    Locus_ID = factor(Locus_ID, levels = locus_levels_rev),
    # Map Score into plotting categories
    Category = factor(Score, levels = c(
      "Definitive", "Strong", "Moderate", "Limited", 
      "Contradictory", "Supportive", "No Known", 
      "Conflict", "No Conflict", "Uncurated"
    ))
  )

# 4. Color Palette specification
# Definitive: Forest Green (#2E7D32)
# Strong: Dark Slate Blue (#1E568B)
# Moderate: Sky Blue (#7BAFD4)
# Limited: Golden Yellow (#FBC02D)
# Contradictory: Dark Red/Crimson (#C62828)
# Supportive: Teal / Cyan (#00897B)
# No Known: Purple / Lavender (#8E24AA)
# Conflict: Dark Red (#C62828)
# No Conflict: Grayish Green / Sage (#A5D6A7)
# Uncurated: Pure White / Light Gray outline (#FFFFFF)

color_palette <- c(
  "Definitive" = "#2E7D32",
  "Strong" = "#1E568B",
  "Moderate" = "#7BAFD4",
  "Limited" = "#FBC02D",
  "Contradictory" = "#C62828",
  "Supportive" = "#00897B",
  "No Known" = "#8E24AA",
  "Conflict" = "#C62828",
  "No Conflict" = "#C8E6C9",
  "Uncurated" = "#FAFAFA"
)

# 5. Build ggplot Heatmap
cat("Generating publication-ready Heatmap Plot...\n")

p <- ggplot(heatmap_df, aes(x = Group_Label, y = Locus_ID, fill = Category)) +
  geom_tile(color = "#E0E0E0", linewidth = 0.35, width = 0.95, height = 0.95) +
  scale_fill_manual(values = color_palette, drop = FALSE, name = "Classification Tier / Conflict Status") +
  scale_x_discrete(position = "top") +
  theme_minimal(base_family = "sans") +
  theme(
    axis.title = element_blank(),
    axis.text.x.top = element_text(face = "bold", size = 9, color = "#222222", vjust = 0.5, lineheight = 1.1),
    axis.text.y = element_text(size = 7.5, color = "#222222", hjust = 1, family = "mono"),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "#CCCCCC", fill = NA, linewidth = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 9.5),
    legend.text = element_text(size = 8.5),
    legend.key.size = unit(0.45, "cm"),
    legend.margin = margin(t = 8, b = 4),
    plot.margin = margin(t = 12, r = 16, b = 12, l = 16)
  ) +
  guides(
    fill = guide_legend(nrow = 1, byrow = TRUE, title.position = "left", title.vjust = 0.8)
  )

# Save PDF and PNG
pdf_file <- "Figure_2_matched_heatmap.pdf"
png_file <- "Figure_2_matched_heatmap.png"

cat("Saving Figure 2 heatmap to PDF and PNG...\n")
ggsave(pdf_file, plot = p, width = 9, height = 12, units = "in")
ggsave(png_file, plot = p, width = 9, height = 12, dpi = 300, units = "in")

cat("Successfully saved:\n")
cat(" -", pdf_file, "\n")
cat(" -", png_file, "\n")
