# 07_generate_supplemental_file_3.R
# Generates Supplemental_File_3_criTRia_Dataset_R2.tsv
# from Supplemental_File_4_criTRia-curations_R1.tsv and matched_gencc_curations.tsv.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

# 1. Load Supplemental File 4
critria_files <- list.files(pattern = "Supplemental_File_4_criTRia-curations.*\\.tsv")
if (length(critria_files) == 0) stop("Supplemental_File_4 not found in current directory.")
critria_file <- critria_files[1]
critria <- read_tsv(critria_file, show_col_types = FALSE)

# 2. Load Matched GenCC Curations
matched_file <- "matched_gencc_curations.tsv"
if (!file.exists(matched_file)) stop("matched_gencc_curations.tsv not found in current directory.")
matched <- read_tsv(matched_file, show_col_types = FALSE)

# Standardize classifications helper
standardize_criTRia_class <- function(x) {
  case_when(
    x %in% c("Disputed", "Refuted", "Contradictory") ~ "Contradictory",
    x == "Definitive" ~ "Definitive",
    x == "Strong" ~ "Strong",
    x == "Moderate" ~ "Moderate",
    x == "Limited" ~ "Limited",
    TRUE ~ x
  )
}

# Target submitters
target_groups <- c("Ambry", "ClinGen", "G2P", "Illumina", "Lab MM", "Labcorp", "Myriad", "Orphanet", "PanelApp")

# A. Add criTRia rows (all 65 loci)
critria_rows <- critria %>%
  transmute(
    Gene = Locus_ID,
    Group = "criTRia",
    categorical_score = standardize_criTRia_class(classification)
  )

# B. Add matched GenCC rows
matched_rows <- matched %>%
  filter(Submitter %in% target_groups) %>%
  transmute(
    Gene = Locus_ID,
    Group = Submitter,
    categorical_score = Score
  )

# Combine and Sort: Primary by Gene (alphabetical), Secondary with criTRia first, then other submitters alphabetically
combined_df <- bind_rows(critria_rows, matched_rows) %>%
  mutate(Group_Rank = ifelse(Group == "criTRia", 0, 1)) %>%
  arrange(Gene, Group_Rank, Group) %>%
  select(Gene, Group, categorical_score)

# 3. Write TSV
tsv_out <- "Supplemental_File_3_criTRia_Dataset_R2.tsv"
write_tsv(combined_df, tsv_out)

cat("Successfully generated:", tsv_out, "(", nrow(combined_df), "rows across", n_distinct(combined_df$Gene), "loci)\n")
