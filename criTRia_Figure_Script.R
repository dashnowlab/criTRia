#set linkgroup= [something]

library(tidyverse)
library(ggplot2)
library(gt)

library(ggplot2) # Need ggplot2 loaded first
library(cowplot)
theme_set(theme_cowplot()) # Sets the default for subsequent plots 

##read and prep data
data <- read_csv(file="criTRia_Dataset.csv")
data %>% rename("categorical_score" = any_of("Categorical Score")) -> data

# Read in the updated STRchive-scored loci
read_tsv("criTRia-curations.tsv", col_select = c("Locus_ID", "Source", "classification")) %>%
  filter(Source == "criTRia") %>%
  rename("Gene" = "Locus_ID", "Group" = "Source", "categorical_score" = "classification") -> strchive_data
# group Refuted and Disputed under Contradictory to simplify plotting
strchive_data %>% mutate(categorical_score = case_when(
  categorical_score == "Refuted" ~ "Contradictory",
  categorical_score == "Disputed" ~ "Contradictory",
  TRUE ~ categorical_score  # Keep all other values as they are
)) -> strchive_data

# Replace all criTRia scores with updated ones
data %>%
  filter(Group != "criTRia") %>%
  bind_rows(strchive_data) -> data

# Set factor levels
data$categorical_score <- factor(data$categorical_score,levels = c("Contradictory","Limited","Moderate","Strong","Supportive","Definitive"),ordered = TRUE)

# Save updated data
write_csv(data, "criTRia_Dataset.csv")


unique(data$categorical_score)
data %>% filter(
  Group %in% c(
    "Clingen",
    "Ambry",
    "Labcorp", 
    "criTRia",
    "G2P"
    )
  ) %>% group_by(
    Gene
    ) %>% summarise(
      num_groups_present = n_distinct (
        Group
        )
      ) %>% filter(
        num_groups_present >= 4
        ) -> shared.genes

##Jitter
jitterPlot <- ggplot(
  data = data,
  aes(
    x=Gene,
    y=Group,
    color=categorical_score,
    size=categorical_score
    )
  ) +  
  geom_jitter(
    width = 0.2,
    height = 0.2
    )
jitterPlot
ggsave("jitter_plot.pdf", plot = jitterPlot, width = 12, height = 8)
##Bar
#criTRia vs GeneCC
plot_data <- data %>%
  count(categorical_score, Group) %>%
  complete(
    categorical_score,
    Group,
    fill = list(n = 0)
  )
criTRiaVsGeneCC<-ggplot(
  data = plot_data,
  aes(
    x = categorical_score,
    y = n,
    fill = Group
    )
  )+
  geom_col(
    position = position_dodge2(
      preserve = "single",
      padding = 0.1
      ),
    width = 0.7
  )+
  geom_vline(
    xintercept = seq(1.5, length(unique(plot_data$categorical_score)) - 0.5, by = 1),
    color = "grey95",
    linewidth = 0.5
  )+
  scale_x_discrete(
    expand = expansion(
      mult = c(
        0.01, 0.01
        )
      )
    )+
  labs(
      title = "GeneCC Scoring Vs. criTRia",
      x = "Categorical Score",
      y = "Number of Genes Scored"
    )+
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9),
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.4, "cm")
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,        # from option 3
      byrow = TRUE,
      title.position = "left"
    )
  )
criTRiaVsGeneCC
ggsave("criTRia_vs_genecc_barplot.pdf", plot = criTRiaVsGeneCC, width = 10, height = 6)

criTRiaPlot <- ggplot(
  data = data,
  aes(
    x=categorical_score
    )
  )+
  geom_bar(
    position = "dodge",
    width = 0.4,
    fill = "#83bbc3"
    )+
  labs(
    title = "criTRia Scoring",
    x = "Categorical Score",
    y = "Number of Genes Scored")
criTRiaPlot
ggsave("criTRia_scoring_barplot.pdf", plot = criTRiaPlot, width = 8, height = 6)

##UpSet Plot
library(UpSetR)
library(data.table)
dt <- fread("criTRia_Dataset.csv")
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
# Drop duplicates of Gene–Group pairs
dt_unique <- unique(dt[, .(gene, group)])
# Build list of genes per group
sets <- split(dt_unique$gene, dt_unique$group)
sets <- lapply(sets, unique)
# Convert to incidence matrix
incidence <- UpSetR::fromList(sets)
# Plot
# Open the file device (e.g., pdf, png, tiff)
pdf("upset_plot.pdf", width = 10, height = 7)
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
                               levels = c("Definitive", "Supportive", "Strong", "Moderate", "Limited", "Contradictory"))
# Order by number of associations scored
data$Group = factor(data$Group,
                    levels = names(sort(table(data$Group), decreasing = T))
)

my_colors <- c(
  "Definitive" = "#59A14F",
  "Supportive" = "#8CD17D",
  "Strong" = "#4E79A7",
  "Moderate" =  "#A0CBE8",
  "Limited" =  "#EDC948",
  "Contradictory" = "#F28E2B"
)

group_counts = data.frame(Group = names(sort(table(data$Group), decreasing = T)), 
                          Count = as.numeric(sort(table(data$Group), decreasing = T)))

counts <- group_counts$Count[match(levels(data$Group), group_counts$Group)]

heatmapPlot <- ggplot(data, aes(y = Gene, x = Group)) +
  geom_tile(aes(fill = categorical_score)) + 
  scale_fill_manual(values = my_colors) +
  scale_x_discrete(sec.axis = dup_axis(labels = counts, name = "Count")) +
  labs(title = "criTRia Vs. GeneCC Scoring", x = "Group", y = "Gene") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"), axis.text = element_text(size = 10))
heatmapPlot
ggsave('heatmap.pdf', plot = heatmapPlot, height = 20, width = 12)
