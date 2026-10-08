# 03_figures.R
# Figures comparing criTRia and GenCC classifications.
# Run from the repository root (see run_all.sh).
# Input:  paper/supp3_dataset.tsv (from scripts/01_download_and_match.py)
# Paper:  paper/fig2_heatmap.pdf, paper/fig2_heatmap.png (Figure 2)
# Not in the paper (exploratory/): jitter_plot.pdf, criTRia_vs_genecc_barplot.pdf,
#   criTRia_scoring_barplot.pdf, upset_plot.pdf, tier_hierarchy.pdf/.png

library(tidyverse)
library(ggplot2)
library(gt)
library(ggnewscale)

library(ggplot2) # Need ggplot2 loaded first
library(cowplot)
theme_set(theme_cowplot()) # Sets the default for subsequent plots 

##read and prep data
data <- read_tsv("paper/supp3_dataset.tsv", show_col_types = FALSE)

# Set factor levels
data$categorical_score <- factor(data$categorical_score,levels = c("Contradictory","Limited","Moderate","Strong","Supportive","Definitive","No Known"),ordered = TRUE)

my_colors <- c(
  "Definitive" = "#59A14F",
  "Strong" = "#4E79A7",
  "Moderate" =  "#A0CBE8",
  "Supportive" = "#D6EAF8",
  "Limited" =  "#EDC948",
  "Contradictory" = "#F28E2B",
  "No Known" = "#D3D3D3"
)

score_sizes <- c(
  "Contradictory" = 2, "Limited" = 2.5, "Moderate" = 3,
  "Strong" = 3.5, "Supportive" = 4, "Definitive" = 4.5,
  "No Known" = 2
)

unique(data$categorical_score)

##Jitter
jitterPlot <- ggplot(
  data = data,
  aes(
    x=Locus_ID,
    y=Group,
    color=categorical_score,
    size=categorical_score
    )
  ) +
  geom_jitter(
    width = 0.2,
    height = 0.2
    ) +
  scale_color_manual(values = my_colors) +
  scale_size_manual(values = score_sizes)
ggsave("exploratory/jitter_plot.pdf", plot = jitterPlot, width = 12, height = 8)
##Bar
#criTRia vs GeneCC
group_totals <- data %>%
  count(Group, name = "total") %>%
  arrange(desc(total))

plot_data <- data %>%
  count(categorical_score, Group) %>%
  complete(categorical_score, Group, fill = list(n = 0)) %>%
  mutate(Group = factor(Group, levels = group_totals$Group))

criTRiaVsGeneCC <- ggplot(
  data = plot_data,
  aes(x = Group, y = n, fill = categorical_score)
  ) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = my_colors, name = "Categorical Score") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "GenCC Scoring Vs. criTRia",
    x = "Scoring Group",
    y = "Number of Loci Scored"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9),
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.4, "cm"),
    axis.text.x = element_text(angle = 35, hjust = 1)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, title.position = "left"))
ggsave("exploratory/criTRia_vs_genecc_barplot.pdf", plot = criTRiaVsGeneCC, width = 10, height = 6)

criTRiaPlot <- ggplot(
  data = subset(data, Group == "criTRia"),
  aes(
    x = categorical_score,
    fill = categorical_score
    )
  )+
  geom_bar(
    position = "dodge",
    width = 0.4
    )+
  scale_fill_manual(values = my_colors, guide = "none") +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 4) +
  labs(
    title = "criTRia Scoring",
    x = "Categorical Score",
    y = "Number of Genes Scored")
ggsave("exploratory/criTRia_scoring_barplot.pdf", plot = criTRiaPlot, width = 8, height = 6)

table(subset(data, Group == "criTRia")$categorical_score)

##UpSet Plot
library(UpSetR)
library(data.table)
dt <- fread("paper/supp3_dataset.tsv")
# Make sure column names are consistent
setnames(
  dt, 
  old = names(
    dt),
  new = tolower(
    names(
      dt
      )
    )
  )
# Drop duplicates of locus–group pairs
dt_unique <- unique(dt[, .(locus_id, group)])
# Build list of loci per group
sets <- split(dt_unique$locus_id, dt_unique$group)
sets <- lapply(sets, unique)
# Convert to incidence matrix
incidence <- UpSetR::fromList(sets)
# Plot
# Open the file device (e.g., pdf, png, tiff)
pdf("exploratory/upset_plot.pdf", width = 10, height = 7)
upset(
  incidence,
  nsets = min(10, ncol(incidence)),   # show up to 10 groups
  nintersects = 19,                   # show top 30 intersections
  order.by = c("degree", "freq"),
  decreasing = c(TRUE, TRUE),
  empty.intersections = "on",
  mainbar.y.label = "Intersection size",
  sets.x.label = "Genes per group"
)

# Close/save plot
dev.off()


##Heat map
# Change order of some values
data$categorical_score = factor(data$categorical_score,
  levels = c("Definitive", "Strong", "Moderate", "Supportive", "Limited", "Contradictory", "No Known"))
# Order by number of associations scored
data$Group = factor(data$Group,
                    levels = names(sort(table(data$Group), decreasing = T))
)

group_counts = data.frame(Group = names(sort(table(data$Group), decreasing = T)),
                          Count = as.numeric(sort(table(data$Group), decreasing = T)))

counts <- group_counts$Count[match(levels(data$Group), group_counts$Group)]

# Sort locus axis by the gene symbol (part after the last underscore)
gene_order <- data %>%
  distinct(Locus_ID) %>%
  mutate(gene_symbol = sub(".*_", "", Locus_ID)) %>%
  arrange(desc(gene_symbol)) %>%
  pull(Locus_ID)
data$Locus_ID <- factor(data$Locus_ID, levels = gene_order)

# Discordance flags (see discordance.R)
source("scripts/discordance.R")
flag_groups <- c("vs criTRia", "Within GenCC")
conflict_data <- discordance_flags(data) %>%
  rename(`vs criTRia` = vs_criTRia, `Within GenCC` = within_GenCC)

n_groups <- nlevels(data$Group)
n_flags <- sapply(flag_groups, function(f) sum(conflict_data[[f]]))
data$Group <- factor(data$Group, levels = c(levels(data$Group), flag_groups))

conflict_col <- conflict_data %>%
  select(Locus_ID, all_of(flag_groups)) %>%
  pivot_longer(all_of(flag_groups), names_to = "Group", values_to = "flagged") %>%
  mutate(
    Locus_ID = factor(Locus_ID, levels = levels(data$Locus_ID)),
    Group = factor(Group, levels = levels(data$Group)),
    conflict = ifelse(flagged, as.character(Group), "No")
  )

heatmapPlot <- ggplot(data, aes(y = Locus_ID, x = Group)) +
  geom_tile(aes(fill = categorical_score)) +
  scale_fill_manual(values = my_colors, name = "Categorical Score") +
  new_scale_fill() +
  geom_tile(
    data = conflict_col,
    aes(y = Locus_ID, x = Group, fill = conflict)
  ) +
  scale_fill_manual(
    values = c("vs criTRia" = "#E15759", "Within GenCC" = "#B07AA1", "No" = "#F5F5F5"),
    name = "Evidence Tier\nConflict"
  ) +
  geom_vline(xintercept = n_groups + 0.5, color = "grey50", linewidth = 0.5, linetype = "dashed") +
  scale_x_discrete(
    labels = function(x) ifelse(x %in% flag_groups, sub(" ", "\n", x), x),
    sec.axis = dup_axis(labels = c(counts, n_flags), name = "Count")
  ) +
  labs(title = "criTRia Vs. GenCC Scoring", x = "Scoring Group", y = "Disease and Gene") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

# Tier hierarchy legend: colored dots with Supportive range box (Limited–Definitive)
tier_legend_tiers <- c("Definitive", "Strong", "Moderate", "Limited", "Contradictory", "No Known")
tier_legend_df <- data.frame(
  category = factor(tier_legend_tiers, levels = rev(tier_legend_tiers)),
  x = 0
)

tierLegend <- ggplot(tier_legend_df, aes(x = x, y = category)) +
  annotate("rect",
           xmin = 0.4, xmax = 1.7,
           ymin = 2.5, ymax = 6.5,
           fill = "#D6EAF8", color = NA) +
  geom_point(aes(color = category), size = 5) +
  scale_color_manual(values = my_colors[tier_legend_tiers], guide = "none") +
  annotate("text",
           x = 0.25, y = 4.5,
           label = "Supportive",
           color = "black", size = 3.5, hjust = -.15) +
  scale_x_continuous(limits = c(-0.5, 2.0)) +
  labs(title = "Categorical Score") +
  theme_minimal() +
  theme(
    plot.title  = element_text(size = 13, face = "bold", hjust = 10),
    axis.title  = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks  = element_blank(),
    axis.text.y = element_text(size = 11, face = "bold"),
    panel.grid  = element_blank()
  )

combinedHeatmap <- ggdraw() +
  draw_plot(heatmapPlot,  x = 0,    y = 0,    width = 0.83, height = 1) +
  draw_plot(tierLegend,   x = 0.83, y = 0.72, width = 0.17, height = 0.15)
ggsave('paper/fig2_heatmap.pdf', plot = combinedHeatmap, height = 20, width = 16)
ggsave('paper/fig2_heatmap.png', plot = combinedHeatmap, height = 20, width = 16)


##Evidence Tier Hierarchy — shows that criTRia's "Supportive" maps broadly
## across GenCC's Limited, Moderate, Strong, and Definitive categories
gencc_tiers <- c("Definitive", "Strong", "Moderate", "Limited", "Contradictory", "No Known")

tier_df <- data.frame(
  category = factor(gencc_tiers, levels = rev(gencc_tiers)),
  x = 0
)

# Supportive overlaps Limited (y=3) through Definitive (y=6)
tierHierarchyPlot <- ggplot(tier_df, aes(x = x, y = category)) +
  annotate("rect",
           xmin = -0.5, xmax = 2.0,
           ymin = 2.5, ymax = 6.5,
           fill = "#D6EAF8", color = NA) +
  geom_point(aes(color = category), size = 10) +
  scale_color_manual(values = my_colors[gencc_tiers], guide = "none") +
  annotate("text",
           x = 0.7, y = 4,
           label = "Supportive\n(Broad Positive Range)",
           color = "#5B9BD5", size = 4.5, hjust = 0) +
  scale_x_continuous(limits = c(-0.8, 2.2)) +
  labs(title = 'Evidence Tier Hierarchy and "Supportive" Range') +
  theme_minimal() +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title   = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 12),
    panel.grid   = element_blank()
  )
ggsave('exploratory/tier_hierarchy.pdf', plot = tierHierarchyPlot, width = 6, height = 5)
ggsave('exploratory/tier_hierarchy.png', plot = tierHierarchyPlot, width = 6, height = 5)
