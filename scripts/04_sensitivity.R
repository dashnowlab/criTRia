# 04_sensitivity.R
# Score perturbation sensitivity analysis across the criTRia-curated loci: how many
# classifications change if every total score is one point lower or higher.
# Run from the repository root (see run_all.sh).
# Input:   paper/supp4_criTRia_curations.tsv
# Outputs: paper/suppX_sensitivity_analysis.tsv (Supplemental File X, one row per locus)
#          data/processed/sensitivity_analysis_summary.tsv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# 1. Load the criTRia curations
critria_file <- "paper/supp4_criTRia_curations.tsv"
critria <- read_tsv(critria_file, show_col_types = FALSE)

cat("Loaded criTRia curations:", nrow(critria), "loci from", critria_file, "\n\n")

# 2. Classification rule, matching summarize_curations() in STRchive's
# scripts/check-curations.py:
# - Disputed and Refuted are set manually and never change.
# - The total score (capped at 18) is rounded to the nearest integer, with halves
#   rounded to even as in Python's round() (R's round() does the same), so 11.5 -> 12.
# - Rounded 12-18: Definitive if >= 2 publications at least 3 years apart, else Strong.
# - Rounded 7-11: Moderate. Rounded 1-6: Limited. Rounded 0: No Known Relationship.
classify_criTRia <- function(score, pub_int, pub_count, orig_class) {
  if (orig_class %in% c("Disputed", "Refuted", "Contradictory")) {
    return(orig_class)
  }
  rounded <- round(min(score, 18))
  pub_int <- ifelse(is.na(pub_int), 0, pub_int)
  if (rounded >= 12) {
    if (pub_count >= 2 && pub_int >= 3) "Definitive" else "Strong"
  } else if (rounded >= 7) {
    "Moderate"
  } else if (rounded >= 1) {
    "Limited"
  } else {
    "No Known Relationship"
  }
}

# Check the rule reproduces every current classification before perturbing
unperturbed <- mapply(
  classify_criTRia,
  critria$total_score,
  critria$publication_interval_years,
  critria$publication_count,
  critria$classification,
  USE.NAMES = FALSE
)
mismatched <- critria$Locus_ID[unperturbed != critria$classification]
if (length(mismatched) > 0) {
  stop("Classification rule does not reproduce current classifications for: ", paste(mismatched, collapse = ", "))
}

# 3. Calculate Perturbed Classifications
# Every total score is shifted down and up by one point, about the median difference between two
# curators (see 05_curator_agreement.R). A locus is sensitive if either shift changes its classification.
deltas <- c(-1, 1)

# Build per-locus data frame
per_locus_df <- critria %>%
  select(
    Locus_ID,
    Gene,
    Disease_ID,
    Original_Score = total_score,
    Publication_Count = publication_count,
    Publication_Interval_Years = publication_interval_years,
    Original_Classification = classification
  )

# Compute new score and new tier for each delta
for (d in deltas) {
  delta_str <- ifelse(d > 0, paste0("+", d), as.character(d))
  scores_d <- pmax(0, per_locus_df$Original_Score + d)
  classes_d <- mapply(
    classify_criTRia,
    scores_d,
    per_locus_df$Publication_Interval_Years,
    per_locus_df$Publication_Count,
    per_locus_df$Original_Classification,
    USE.NAMES = FALSE
  )
  per_locus_df[[paste0("Score_", delta_str)]] <- scores_d
  per_locus_df[[paste0("Class_", delta_str)]] <- classes_d
}
per_locus_df <- per_locus_df %>%
  mutate(Sensitive = `Class_-1` != Original_Classification | `Class_+1` != Original_Classification)

# 4. Generate Summary Table
# One row each for -1, +1 and either (+/-1): number of loci unchanged and changed, with the changes listed
describe_shifts <- function(to) {
  from <- per_locus_df$Original_Classification
  changed <- from != to
  if (!any(changed)) return("None")
  tibble(shift = paste0(from, " -> ", to), locus = per_locus_df$Locus_ID)[changed, ] %>%
    group_by(shift) %>%
    summarise(n = n(), loci = paste(locus, collapse = ", "), .groups = "drop") %>%
    mutate(desc = paste0(shift, " (n=", n, ": ", loci, ")")) %>%
    pull(desc) %>%
    paste(collapse = "; ")
}
summarise_shift <- function(label, changed, details) {
  tibble(
    Score_Perturbation = label,
    Total_Loci = length(changed),
    Stable_Loci = sum(!changed),
    Shifted_Loci = sum(changed),
    Stability_Percentage = paste0(round(mean(!changed) * 100, 1), "%"),
    Details_of_Reclassified_Loci = details
  )
}
summary_df <- bind_rows(
  summarise_shift("-1 pt", per_locus_df$`Class_-1` != per_locus_df$Original_Classification,
                  describe_shifts(per_locus_df$`Class_-1`)),
  summarise_shift("+1 pt", per_locus_df$`Class_+1` != per_locus_df$Original_Classification,
                  describe_shifts(per_locus_df$`Class_+1`)),
  summarise_shift("+/-1 pt (either)", per_locus_df$Sensitive,
                  paste(per_locus_df$Locus_ID[per_locus_df$Sensitive], collapse = ", "))
)

# 5. Output TSV Files
per_tsv <- "paper/suppX_sensitivity_analysis.tsv"
sum_tsv <- "data/processed/sensitivity_analysis_summary.tsv"

write_tsv(per_locus_df, per_tsv)
write_tsv(summary_df, sum_tsv)

cat("Generated Sensitivity Analysis TSV deliverables:\n")
cat(" -", per_tsv, "\n")
cat(" -", sum_tsv, "\n\n")

# Display Summary in Console
cat("=== Score Perturbation Sensitivity Analysis Summary ===\n\n")
print(as.data.frame(summary_df))
cat("\n")
