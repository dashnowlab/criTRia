# Inter-rater agreement between two criTRia curators, from the per-locus table written by
# scripts/compare_curators.py.
#   - ICC(1,1): one-way random effects, single rater (Shrout & Fleiss 1979), on the calculated
#     scores (Manual_evidence_level does not change scores). One-way because the "group" rater is
#     a different curator from locus to locus. ICC(A,1) (two-way, absolute agreement) is also
#     reported as a sensitivity analysis. Analytic 95% CIs, as in irr::icc.
#   - Weighted Cohen's kappa on ordered classifications: quadratic weights (primary) and linear.
#     95% CI from a percentile bootstrap over loci.
#   - Differences in total score: mean difference (t-test CI), Bland-Altman limits of agreement,
#     median and maximum absolute difference, and proportion of loci within 1.5 points.
#
# Usage: Rscript scripts/05_curator_agreement.R [comparison.tsv] [rater_a] [rater_b] [out.tsv] [figure_prefix] [label_a] [label_b]
# The figure (A: total score per locus; B: agreement estimates) is written to <figure_prefix>.pdf/png.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
in_file  <- if (length(args) >= 1) args[1] else "data/curators/curator-comparison.tsv"
rater_a  <- if (length(args) >= 2) args[2] else "gabriel"
rater_b  <- if (length(args) >= 3) args[3] else "group"
out_file <- if (length(args) >= 4) args[4] else "data/curators/curator-agreement.tsv"

n_boot <- 10000
set.seed(20261027)

data <- read_tsv(in_file, show_col_types = FALSE)
col <- function(rater, field) data[[paste0(rater, "_", field)]]

# Leave out loci that either curator has not finished: an evidence entry with a source or details
# but no score (flagged by compare_curators.py), or a missing score/classification. They are
# included automatically once the curation is filled in and the comparison table is rebuilt.
needed <- c("incomplete_entries", "total_score", "genetic_score", "experimental_score",
            "classification_calculated", "classification_final")
missing <- setdiff(paste0(rep(c(rater_a, rater_b), each = length(needed)), "_", needed), names(data))
if (length(missing) > 0) stop("Columns missing from ", in_file, ": ", paste(missing, collapse = ", "))
unfinished <- rep(FALSE, nrow(data))
for (rater in c(rater_a, rater_b)) {
  incomplete <- col(rater, "incomplete_entries")
  unfinished_rater <- (!is.na(incomplete) & incomplete != "") |
    rowSums(is.na(data[paste0(rater, "_", needed[-1])])) > 0
  for (i in which(unfinished_rater)) {
    message("Excluding ", data$Locus_ID[i], ": ", rater, " unfinished",
            if (!is.na(incomplete[i]) && incomplete[i] != "") paste0(" (no score for ", incomplete[i], ")") else " (missing score or classification)")
  }
  unfinished <- unfinished | unfinished_rater
}
data <- data[!unfinished, ]
message(nrow(data), " loci completed by both ", rater_a, " and ", rater_b)

# Classification scales, lowest to highest.
# Primary: as in 03_figures.R, Disputed and Refuted are merged into Contradictory, above No Known.
# Alternative: Refuted and Disputed kept separate, below Limited.
merge_contradictory <- function(x) ifelse(x %in% c("Disputed", "Refuted"), "Contradictory", x)
scales <- list(
  final_contradictory = list(
    field = "classification_final", recode = merge_contradictory,
    levels = c("No Known Relationship", "Contradictory", "Limited", "Moderate", "Strong", "Definitive")),
  final_disputed_refuted_separate = list(
    field = "classification_final", recode = identity,
    levels = c("No Known Relationship", "Refuted", "Disputed", "Limited", "Moderate", "Strong", "Definitive")),
  calculated = list(
    field = "classification_calculated", recode = identity,
    levels = c("No Known Relationship", "Limited", "Moderate", "Strong", "Definitive"))
)

# ICC(1,1) with F-test confidence interval (same formulas as irr::icc, model = "oneway")
icc_1 <- function(x, y, alpha = 0.05) {
  ratings <- cbind(x, y)
  n <- nrow(ratings); k <- ncol(ratings)
  grand <- mean(ratings)
  msb <- k * sum((rowMeans(ratings) - grand)^2) / (n - 1)
  msw <- sum((ratings - rowMeans(ratings))^2) / (n * (k - 1))
  icc <- (msb - msw) / (msb + (k - 1) * msw)
  f0 <- msb / msw
  fl <- f0 / qf(1 - alpha / 2, n - 1, n * (k - 1))
  fu <- f0 * qf(1 - alpha / 2, n * (k - 1), n - 1)
  c(estimate = icc, lower = (fl - 1) / (fl + k - 1), upper = (fu - 1) / (fu + k - 1))
}

# ICC(A,1) with McGraw & Wong confidence interval (same formulas as irr::icc)
icc_a1 <- function(x, y, alpha = 0.05) {
  ratings <- cbind(x, y)
  n <- nrow(ratings); k <- ncol(ratings)
  grand <- mean(ratings)
  ss_rows <- k * sum((rowMeans(ratings) - grand)^2)
  ss_cols <- n * sum((colMeans(ratings) - grand)^2)
  ss_err  <- sum((ratings - grand)^2) - ss_rows - ss_cols
  msr <- ss_rows / (n - 1)
  msc <- ss_cols / (k - 1)
  mse <- ss_err / ((n - 1) * (k - 1))
  icc <- (msr - mse) / (msr + (k - 1) * mse + k / n * (msc - mse))
  a <- k * icc / (n * (1 - icc))
  b <- 1 + k * icc * (n - 1) / (n * (1 - icc))
  v <- (a * msc + b * mse)^2 / ((a * msc)^2 / (k - 1) + (b * mse)^2 / ((n - 1) * (k - 1)))
  fl <- qf(1 - alpha / 2, n - 1, v)
  fu <- qf(1 - alpha / 2, v, n - 1)
  lower <- n * (msr - fl * mse) / (fl * (k * msc + (k * n - k - n) * mse) + n * msr)
  upper <- n * (fu * msr - mse) / (k * msc + (k * n - k - n) * mse + n * fu * msr)
  c(estimate = icc, lower = lower, upper = upper)
}

# Cohen's weighted kappa for two raters on an ordered scale
weighted_kappa <- function(x, y, levels, weights = c("linear", "quadratic")) {
  weights <- match.arg(weights)
  m <- length(levels)
  observed <- table(factor(x, levels = levels), factor(y, levels = levels)) / length(x)
  expected <- outer(rowSums(observed), colSums(observed))
  distance <- abs(outer(seq_len(m), seq_len(m), "-")) / (m - 1)
  disagreement <- if (weights == "linear") distance else distance^2
  1 - sum(disagreement * observed) / sum(disagreement * expected)
}

boot_ci <- function(stat, n, alpha = 0.05) {
  estimates <- replicate(n_boot, {
    idx <- sample.int(n, replace = TRUE)
    stat(idx)
  })
  quantile(estimates, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE, names = FALSE)
}

results <- list()
add_result <- function(...) results[[length(results) + 1]] <<- tibble(...)

# ICC on scores: all loci, then excluding loci with a manual classification override
overridden <- !is.na(col(rater_a, "manual_evidence_level")) | !is.na(col(rater_b, "manual_evidence_level"))
for (subset in c("all", "no_manual_override")) {
  keep <- if (subset == "all") rep(TRUE, nrow(data)) else !overridden
  for (score in c("total_score", "genetic_score", "experimental_score")) {
    for (model in c("ICC(1,1)", "ICC(A,1)")) {
      icc_fun <- if (model == "ICC(1,1)") icc_1 else icc_a1
      res <- icc_fun(col(rater_a, score)[keep], col(rater_b, score)[keep])
      add_result(statistic = model, measure = score, subset = subset, n = sum(keep),
                 estimate = res[["estimate"]], ci_lower = res[["lower"]], ci_upper = res[["upper"]],
                 ci_method = "F-test")
    }
  }

  # Differences in total score (rater_a minus rater_b): bias with t-based CI, Bland-Altman limits of
  # agreement, and the size of typical and largest differences
  diff <- col(rater_a, "total_score")[keep] - col(rater_b, "total_score")[keep]
  ci <- t.test(diff)$conf.int
  add_result(statistic = "Mean difference", measure = "total_score", subset = subset, n = length(diff),
             estimate = mean(diff), ci_lower = ci[1], ci_upper = ci[2], ci_method = "t-test")
  add_result(statistic = "95% limits of agreement", measure = "total_score", subset = subset, n = length(diff),
             estimate = NA, ci_lower = mean(diff) - 1.96 * sd(diff), ci_upper = mean(diff) + 1.96 * sd(diff),
             ci_method = "Bland-Altman (mean ± 1.96 SD)")
  add_result(statistic = "Median absolute difference", measure = "total_score", subset = subset, n = length(diff),
             estimate = median(abs(diff)), ci_lower = NA, ci_upper = NA, ci_method = NA)
  add_result(statistic = "Maximum absolute difference", measure = "total_score", subset = subset, n = length(diff),
             estimate = max(abs(diff)), ci_lower = NA, ci_upper = NA, ci_method = NA)
  add_result(statistic = "Within 1.5 points (proportion)", measure = "total_score", subset = subset, n = length(diff),
             estimate = mean(abs(diff) <= 1.5), ci_lower = NA, ci_upper = NA, ci_method = NA)
  largest <- which(abs(diff) == max(abs(diff)))
  message("Largest total score difference (", subset, "): ",
          paste0(data$Locus_ID[keep][largest], " (", diff[largest], ")", collapse = ", "))
}

# Weighted kappa and percent agreement on classifications
for (scale_name in names(scales)) {
  scale <- scales[[scale_name]]
  x <- scale$recode(col(rater_a, scale$field))
  y <- scale$recode(col(rater_b, scale$field))
  unknown <- setdiff(c(x, y), scale$levels)
  if (length(unknown) > 0) stop("Classifications not on the ", scale_name, " scale: ", paste(unknown, collapse = ", "))

  for (w in c("linear", "quadratic")) {
    ci <- boot_ci(function(idx) weighted_kappa(x[idx], y[idx], scale$levels, w), length(x))
    add_result(statistic = paste0("Weighted kappa (", w, ")"), measure = scale_name, subset = "all", n = length(x),
               estimate = weighted_kappa(x, y, scale$levels, w), ci_lower = ci[1], ci_upper = ci[2],
               ci_method = paste0("Percentile bootstrap (", n_boot, ")"))
  }
  ci <- boot_ci(function(idx) mean(x[idx] == y[idx]), length(x))
  add_result(statistic = "Exact agreement", measure = scale_name, subset = "all", n = length(x),
             estimate = mean(x == y), ci_lower = ci[1], ci_upper = ci[2],
             ci_method = paste0("Percentile bootstrap (", n_boot, ")"))

  cat("\nClassification table,", scale_name, "(rows:", rater_a, "; columns:", rater_b, ")\n")
  present <- scale$levels[scale$levels %in% c(x, y)]
  print(table(factor(x, levels = present), factor(y, levels = present), dnn = NULL))
}

results <- bind_rows(results) %>%
  mutate(across(c(estimate, ci_lower, ci_upper), ~ round(.x, 3)))
cat("\n")
print(as.data.frame(results), row.names = FALSE)
write_tsv(results, out_file)
cat("\nWrote", out_file, "\n")

# ---- Figures ----------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})
fig_prefix <- if (length(args) >= 5) args[5] else "paper/figX_curator_agreement"
# How each rater is described in the figure (roles rather than names)
label_a <- if (length(args) >= 6) args[6] else "Independent curator"
label_b <- if (length(args) >= 7) args[7] else "Original curation"

primary_color <- "#2a78d6"
secondary_color <- "#898781"
text_color <- "#0b0b0b"
muted_text <- "#52514e"
grid_color <- "#e1e0d9"

theme_agreement <- theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color = grid_color, linewidth = 0.3),
        axis.text = element_text(color = muted_text),
        axis.title = element_text(color = muted_text),
        strip.text = element_text(color = text_color, face = "bold", hjust = 0),
        plot.title = element_text(color = text_color, face = "bold"),
        plot.subtitle = element_text(color = muted_text))

# Forest plot of the primary analyses. The genetic/experimental breakdown and the sensitivity analyses
# (primary = FALSE) are in the results TSV; set show_sensitivity to TRUE to add them to the plot in grey.
show_sensitivity <- FALSE
forest_rows <- tribble(
  ~statistic,                    ~measure,                          ~subset,              ~panel,                     ~label,                                      ~primary,
  "ICC(1,1)",                    "total_score",                     "all",                "Scores",              "Total score\nICC(1,1)",                           TRUE,
  "ICC(1,1)",                    "total_score",                     "no_manual_override", "Scores",              "  excluding manual overrides",              FALSE,
  "ICC(A,1)",                    "total_score",                     "all",                "Scores",              "  two-way, absolute agreement",             FALSE,
  "ICC(1,1)",                    "genetic_score",                   "all",                "Scores",              "Genetic evidence score",                    FALSE,
  "ICC(1,1)",                    "genetic_score",                   "no_manual_override", "Scores",              "  excluding manual overrides",              FALSE,
  "ICC(A,1)",                    "genetic_score",                   "all",                "Scores",              "  two-way, absolute agreement",             FALSE,
  "ICC(1,1)",                    "experimental_score",              "all",                "Scores",              "Experimental evidence score",               FALSE,
  "ICC(1,1)",                    "experimental_score",              "no_manual_override", "Scores",              "  excluding manual overrides",              FALSE,
  "ICC(A,1)",                    "experimental_score",              "all",                "Scores",              "  two-way, absolute agreement",             FALSE,
  "Weighted kappa (quadratic)",  "final_contradictory",             "all",                "Classification",           "Weighted kappa\n(quadratic)",                  TRUE,
  "Weighted kappa (linear)",     "final_contradictory",             "all",                "Classification",           "  linear weights",                          FALSE,
  "Weighted kappa (quadratic)",  "final_disputed_refuted_separate", "all",                "Classification",           "  Disputed and Refuted separate",           FALSE,
  "Weighted kappa (quadratic)",  "calculated",                      "all",                "Classification",           "  calculated (before manual override)",     FALSE,
  "Exact agreement",             "final_contradictory",             "all",                "Classification",           "Identical\nclassification\n(proportion)",     TRUE
)
forest <- forest_rows %>%
  filter(primary | show_sensitivity) %>%
  left_join(results, by = c("statistic", "measure", "subset")) %>%
  mutate(value_label = sprintf("%.2f (%.2f\u2013%.2f)", estimate, ci_lower, ci_upper),
         analysis = factor(ifelse(primary, "Primary", "Sensitivity"), levels = c("Primary", "Sensitivity")))

# Panel B: one row per estimate, rows positioned by number (so repeated sub-labels stay separate),
# grouped into sections with a bold heading and a hairline between them
section_gap <- 1.3
forest <- forest %>%
  mutate(section = match(panel, unique(forest_rows$panel)),
         y = -(row_number() - 1) - (section - 1) * section_gap)
section_headings <- forest %>% group_by(section, panel) %>% summarise(y = max(y) + 0.75, .groups = "drop")
section_breaks <- section_headings %>% filter(section > 1) %>% mutate(y = y + 0.4)

forest_plot <- ggplot(forest, aes(x = estimate, y = y, color = analysis)) +
  geom_hline(data = section_breaks, aes(yintercept = y), color = grid_color, linewidth = 0.4) +
  geom_text(data = section_headings, aes(x = 0, y = y, label = panel), inherit.aes = FALSE,
            hjust = 0, size = 3.4, fontface = "bold", color = text_color) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0, linewidth = 0.6, orientation = "y") +
  geom_point(aes(size = analysis)) +
  geom_text(aes(x = 1.04, label = value_label), hjust = 0, size = 3.2,
            color = ifelse(forest$primary, text_color, muted_text)) +
  scale_y_continuous(breaks = forest$y, labels = forest$label, expand = expansion(add = c(0.6, 0.3))) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.01))) +
  scale_color_manual(values = c(Primary = primary_color, Sensitivity = secondary_color), name = NULL, drop = FALSE) +
  scale_size_manual(values = c(Primary = 2.8, Sensitivity = 2), name = NULL, drop = FALSE) +
  coord_cartesian(xlim = c(0, 1), clip = "off") +
  labs(x = "Estimate (95% CI)", y = NULL) +
  theme_agreement +
  theme(axis.text.y = element_text(color = text_color, hjust = 0),
        legend.position = if (show_sensitivity) "bottom" else "none")

# Panel A: total score per locus. Dashed line = identical scores; dotted lines = approximate score cutoffs
# between Limited, Moderate and Strong/Definitive (scores are rounded before classifying, so the cutoffs
# fall at 6.5 and 11.5; Strong vs Definitive depends on replication, not score). Points are coloured by
# whether the two final classifications agree, and every disagreement is labelled.
disagree_color <- "#eb6834"
score_limit <- 18
category_cutoffs <- c(6.5, 11.5)
scores <- tibble(Locus_ID = data$Locus_ID, a = col(rater_a, "total_score"), b = col(rater_b, "total_score"),
                 agreement = factor(ifelse(col(rater_a, "classification_final") == col(rater_b, "classification_final"),
                                           "Same classification", "Different classification"),
                                    levels = c("Same classification", "Different classification"))) %>%
  # Unlabelled points get "" rather than NA so the labels are still pushed clear of them
  mutate(label = ifelse(agreement == "Different classification", Locus_ID, "")) %>%
  # Loci with identical scores from both curators would hide each other, so spread them slightly along the
  # diagonal (x and y together), which keeps each point's distance from the line of identical scores
  group_by(a, b) %>%
  mutate(offset = (row_number() - (n() + 1) / 2) * 0.3, a_plot = a + offset, b_plot = b + offset) %>%
  ungroup()
bands <- tibble(mid = c(category_cutoffs[1] / 2, mean(category_cutoffs), (category_cutoffs[2] + score_limit) / 2),
                name = c("Limited", "Moderate", "Strong/\nDefinitive"))

scatter_plot <- ggplot(scores, aes(x = a_plot, y = b_plot)) +
  geom_vline(xintercept = category_cutoffs, color = secondary_color, linewidth = 0.3, linetype = "dotted") +
  geom_hline(yintercept = category_cutoffs, color = secondary_color, linewidth = 0.3, linetype = "dotted") +
  geom_text(data = bands, aes(x = mid, y = score_limit + 1.2, label = name), inherit.aes = FALSE,
            size = 2.7, color = muted_text, lineheight = 0.85, vjust = 0) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = secondary_color, linewidth = 0.4) +
  geom_point(aes(color = agreement, shape = agreement), size = 2.4, alpha = 0.85) +
  geom_text_repel(aes(label = label), size = 2.8, color = muted_text, min.segment.length = 0,
                  box.padding = 0.5, point.padding = 0.4, point.size = 2.4, segment.color = secondary_color,
                  segment.size = 0.3, max.overlaps = Inf, seed = 1) +
  scale_color_manual(values = c("Same classification" = primary_color, "Different classification" = disagree_color),
                     labels = c("Same", "Different"), name = "Classification") +
  scale_shape_manual(values = c("Same classification" = 16, "Different classification" = 17),
                     labels = c("Same", "Different"), name = "Classification") +
  scale_x_continuous(breaks = seq(0, score_limit, 3)) +
  scale_y_continuous(breaks = seq(0, score_limit, 3)) +
  coord_fixed(xlim = c(0, score_limit), ylim = c(0, score_limit), clip = "off") +
  labs(x = paste0(label_a, " total score"), y = paste0(label_b, " total score")) +
  theme_agreement +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.line = element_line(color = grid_color, linewidth = 0.4),
        legend.position = "inside", legend.position.inside = c(1, 0.02), legend.justification = c(1, 0),
        legend.key.height = unit(0.9, "lines"), legend.key.width = unit(0.9, "lines"),
        legend.title = element_text(size = 9, color = text_color),
        legend.text = element_text(size = 9, color = text_color),
        legend.margin = margin(4, 6, 4, 6),
        legend.background = element_rect(fill = "white", color = secondary_color, linewidth = 0.3),
        plot.margin = margin(22, 5.5, 5.5, 5.5))

# Panel B has only a few rows, so it is padded to sit vertically centred beside A rather than stretched
# to its height. It is wrapped so its tag stays level with A's.
figure_height <- 4.7
forest_share <- min(0.9, 0.35 + 0.1 * nrow(forest))
forest_padding <- (1 - forest_share) / 2 * figure_height * 72  # points
figure <- (scatter_plot |
           wrap_elements(full = forest_plot + theme(plot.margin = margin(forest_padding, 95, forest_padding, 10)))) +
  plot_layout(widths = c(1, 1.05)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", color = text_color))
ggsave(paste0(fig_prefix, ".pdf"), figure, width = 10, height = figure_height, device = cairo_pdf)
ggsave(paste0(fig_prefix, ".png"), figure, width = 10, height = figure_height, dpi = 300)
cat("Wrote", paste0(fig_prefix, ".pdf/.png"), "\n")
