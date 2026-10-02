# 08_run_sensitivity_analysis.R
# Executes Score Perturbation Sensitivity Analysis across All 65 criTRia Loci.
# Generates TSV outputs only.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# 1. Locate Supplemental File 4 (Input File) in current directory
critria_files <- list.files(pattern = "Supplemental_File_4_criTRia-curations.*\\.tsv")
if (length(critria_files) == 0) {
  stop("Supplemental_File_4_criTRia-curations file not found in current directory.")
}
critria_file <- critria_files[1]
critria <- read_tsv(critria_file, show_col_types = FALSE)

cat("Loaded criTRia curations:", nrow(critria), "loci from", critria_file, "\n\n")

# 2. Mathematical Classification Logic (criTRia SOP v1)
# For any perturbed score (score + delta, bounded at minimum 0):
# 1. Rule-Based Overrides: If classification is "Disputed" or "Refuted", return original classification.
# 2. Limited: score < 7.0
# 3. Moderate: 7.0 <= score < 12.0
# 4. Strong / Definitive: score >= 12.0
#    - If !is.na(publication_interval_years) and publication_interval_years >= 3.0: "Definitive"
#    - Otherwise: "Strong"

classify_criTRia <- function(score, pub_int, orig_class) {
  if (orig_class %in% c("Disputed", "Refuted", "Contradictory")) {
    return(orig_class)
  }
  if (score < 7.0) {
    return("Limited")
  } else if (score < 12.0) {
    return("Moderate")
  } else {
    # score >= 12.0
    if (!is.na(pub_int) && pub_int >= 3.0) {
      return("Definitive")
    } else {
      return("Strong")
    }
  }
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
sum_tsv <- "sensitivity_analysis_summary.tsv"
per_tsv <- "sensitivity_analysis_per_locus.tsv"

write_tsv(summary_df, sum_tsv)
write_tsv(per_locus_df, per_tsv)

cat("Generated Sensitivity Analysis TSV deliverables:\n")
cat(" -", sum_tsv, "\n")
cat(" -", per_tsv, "\n\n")

# Display Summary in Console
cat("=== Score Perturbation Sensitivity Analysis Summary ===\n\n")
print(as.data.frame(summary_df))
cat("\n")
