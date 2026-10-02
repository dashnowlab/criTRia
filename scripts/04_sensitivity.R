# 04_sensitivity.R
# Score perturbation sensitivity analysis across the criTRia-curated loci.
# Not used in the paper. Run from the repository root (see run_all.sh).
# Input:   paper/supp4_criTRia_curations.tsv
# Outputs: exploratory/sensitivity_analysis_summary.tsv
#          exploratory/sensitivity_analysis_per_locus.tsv

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
deltas <- c(-1.0, -0.5, 0.5, 1.0)

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
  score_col <- paste0("Score_", delta_str)
  class_col <- paste0("Class_", delta_str)
  shift_col <- paste0("Shift_", delta_str)
  
  scores_d <- pmax(0, per_locus_df$Original_Score + d)
  classes_d <- mapply(
    classify_criTRia,
    scores_d,
    per_locus_df$Publication_Interval_Years,
    per_locus_df$Publication_Count,
    per_locus_df$Original_Classification,
    USE.NAMES = FALSE
  )
  
  shifts_d <- ifelse(classes_d == per_locus_df$Original_Classification, "No Shift", paste0(per_locus_df$Original_Classification, " -> ", classes_d))
  
  per_locus_df[[score_col]] <- scores_d
  per_locus_df[[class_col]] <- classes_d
  per_locus_df[[shift_col]] <- shifts_d
}

# 4. Generate Summary Table
# Grouped by Delta, count total loci, loci shifted, stable count, stability percentage, and specific tier shifts
summary_list <- list()

for (d in deltas) {
  delta_str <- ifelse(d > 0, paste0("+", d), as.character(d))
  class_col <- paste0("Class_", delta_str)
  shift_col <- paste0("Shift_", delta_str)
  
  total_loci <- nrow(per_locus_df)
  shifts <- per_locus_df %>% filter(.data[[shift_col]] != "No Shift")
  n_shifted <- nrow(shifts)
  n_stable <- total_loci - n_shifted
  pct_stable <- round((n_stable / total_loci) * 100, 1)
  
  # Format specific shifts summary
  if (n_shifted > 0) {
    shift_details <- shifts %>%
      group_by(.data[[shift_col]]) %>%
      summarise(n = n(), loci = paste(Locus_ID, collapse = ", "), .groups = "drop") %>%
      mutate(desc = paste0(.data[[shift_col]], " (n=", n, ": ", loci, ")")) %>%
      pull(desc) %>%
      paste(collapse = "; ")
  } else {
    shift_details <- "None (100% concordance)"
  }
  
  summary_list[[length(summary_list) + 1]] <- tibble(
    Score_Perturbation = paste0(delta_str, " pt"),
    Total_Loci = total_loci,
    Stable_Loci = n_stable,
    Shifted_Loci = n_shifted,
    Stability_Percentage = paste0(pct_stable, "%"),
    Details_of_Reclassified_Loci = shift_details
  )
}

summary_df <- bind_rows(summary_list)

# 5. Output TSV Files to current directory
sum_tsv <- "exploratory/sensitivity_analysis_summary.tsv"
per_tsv <- "exploratory/sensitivity_analysis_per_locus.tsv"

write_tsv(summary_df, sum_tsv)
write_tsv(per_locus_df, per_tsv)

cat("Generated Sensitivity Analysis TSV deliverables:\n")
cat(" -", sum_tsv, "\n")
cat(" -", per_tsv, "\n\n")

# Display Summary in Console
cat("=== Score Perturbation Sensitivity Analysis Summary ===\n\n")
print(as.data.frame(summary_df))
cat("\n")
