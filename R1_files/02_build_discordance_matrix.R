# analysis/02_build_discordance_matrix.R
# Combines matched GenCC data with criTRia scores into a wide matrix, computes conflicts, and exports Supplemental Table S3b.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

dir.create("analysis/output", recursive = TRUE, showWarnings = FALSE)

# 1. Load data
critria_file <- list.files(pattern = "Supplemental_File_4_criTRia-curations.*\\.tsv")[1]
critria <- read_tsv(critria_file, show_col_types = FALSE)
matched_gencc <- read_tsv("matched_gencc_curations.tsv", show_col_types = FALSE)

# Standardize classifications
standardize_class <- function(x) {
  case_when(
    x == "Definitive" ~ "Definitive",
    x == "Strong" ~ "Strong",
    x == "Moderate" ~ "Moderate",
    x == "Limited" ~ "Limited",
    x %in% c("Disputed Evidence", "Refuted Evidence", "Disputed", "Refuted", "Contradictory") ~ "Contradictory",
    x == "Supportive" ~ "Supportive",
    x %in% c("No Known Disease Relationship", "No Known") ~ "No Known",
    TRUE ~ "Uncurated"
  )
}

# Submitter column order
target_submitters <- c("criTRia", "Labcorp", "Orphanet", "Ambry", "PanelApp", "G2P", "ClinGen", "Lab MM", "Illumina", "Myriad")

# Prepare criTRia entries
critria_clean <- critria %>%
  select(Locus_ID, Gene, Disease_ID, classification, total_score, genetic_evidence_score, experimental_evidence_score, Description) %>%
  mutate(
    Submitter = "criTRia",
    Score = standardize_class(classification)
  )

# Filter matched GenCC to target submitters
matched_filtered <- matched_gencc %>%
  filter(Submitter %in% target_submitters)

# Combine criTRia and GenCC in long format
combined_long <- bind_rows(
  critria_clean %>% select(Locus_ID, Gene, Disease_ID, Submitter, Score),
  matched_filtered %>% select(Locus_ID, Gene, Disease_ID, Submitter, Score)
)

# Pivot wide matrix
matrix_wide <- combined_long %>%
  pivot_wider(id_cols = c(Locus_ID, Gene, Disease_ID), names_from = Submitter, values_from = Score)

# Fill uncurated cells
for (col in target_submitters) {
  if (!col %in% colnames(matrix_wide)) {
    matrix_wide[[col]] <- "Uncurated"
  } else {
    matrix_wide[[col]] <- ifelse(is.na(matrix_wide[[col]]), "Uncurated", matrix_wide[[col]])
  }
}

# Reorder columns
matrix_wide <- matrix_wide %>%
  select(Locus_ID, Gene, Disease_ID, all_of(target_submitters))

# Calculate Conflict:
# "A locus is flagged as Conflict (Red) if at least one group assigned a positive/high category (Definitive, Strong, or Moderate) 
# AND at least one other group assigned a low/negative category (Limited, Contradictory, or No Known).
# Note: The Supportive category (used by Orphanet) is NOT considered when determining conflict.
# If a cell has no matched curation, it is marked as Uncurated and does NOT trigger a conflict."

calc_conflict <- function(row_vals) {
  pos_cats <- c("Definitive", "Strong", "Moderate")
  neg_cats <- c("Limited", "Contradictory", "No Known")
  
  has_pos <- any(row_vals %in% pos_cats)
  has_neg <- any(row_vals %in% neg_cats)
  
  return(has_pos && has_neg)
}

matrix_wide$Conflict <- apply(matrix_wide[, target_submitters], 1, calc_conflict)

# Save wide discordance matrix for downstream use
saveRDS(matrix_wide, "discordance_matrix_wide.rds")
write_tsv(matrix_wide, "discordance_matrix_wide.tsv")

cat("Discordance matrix computed:\n")
cat(" - Total loci:", nrow(matrix_wide), "\n")
cat(" - Discordant (Conflict == TRUE) loci:", sum(matrix_wide$Conflict), "\n")

# Build Supplemental Table S3b for discordant loci
# Required columns: Gene, Locus_ID, Disease_Phenotype, GenCC_Scores, criTRia_Classification, Total_Score, Genetic_Score, Experimental_Score, Discordance_Scenario, Rationale

conflict_loci_df <- matrix_wide %>% filter(Conflict == TRUE)

s3b_list <- list()

for (i in 1:nrow(conflict_loci_df)) {
  loc <- conflict_loci_df$Locus_ID[i]
  g <- conflict_loci_df$Gene[i]
  d_id <- conflict_loci_df$Disease_ID[i]
  
  # criTRia info
  cr_info <- critria %>% filter(Locus_ID == loc)
  cr_class <- standardize_class(cr_info$classification[1])
  tot_sc <- cr_info$total_score[1]
  gen_sc <- cr_info$genetic_evidence_score[1]
  exp_sc <- cr_info$experimental_evidence_score[1]
  desc <- cr_info$Description[1]
  
  # Matched GenCC submitters for this locus
  subs_matched <- matched_filtered %>% filter(Locus_ID == loc)
  
  # Format GenCC_Scores string, e.g. "Labcorp: Strong; PanelApp: Limited; Orphanet: Supportive"
  if (nrow(subs_matched) > 0) {
    gencc_scores_str <- paste(paste0(subs_matched$Submitter, ": ", subs_matched$Score), collapse = "; ")
  } else {
    gencc_scores_str <- "None"
  }
  
  # Determine Discordance Scenario & Rationale
  pos_cats <- c("Definitive", "Strong", "Moderate")
  neg_cats <- c("Limited", "Contradictory", "No Known")
  
  gencc_evals <- subs_matched %>% filter(Score != "Supportive")
  gencc_has_pos <- any(gencc_evals$Score %in% pos_cats)
  gencc_has_neg <- any(gencc_evals$Score %in% neg_cats)
  
  critria_is_pos <- cr_class %in% pos_cats
  critria_is_neg <- cr_class %in% neg_cats
  
  scenario <- ""
  rationale <- ""
  
  if (loc == "DMD_DMD") {
    scenario <- "Locus vs. Gene Discordance (LoF vs. Non-pathogenic Repeat)"
    rationale <- "DMD loss-of-function SNVs/deletions are Definitive causes of Duchenne/Becker muscular dystrophy in GenCC, but the specific intronic GAA repeat expansion is Refuted / Contradictory due to high allele frequency in healthy population controls."
  } else if (loc == "CPEO_POLG") {
    scenario <- "Locus vs. Gene Discordance (Coding SNV vs. Disputed Repeat)"
    rationale <- "POLG point mutations and deletions are Definitive causes of progressive external ophthalmoplegia (CPEO), but evidence for the noncoding/CAG repeat causing CPEO is disputed/lacking locus-specific causality."
  } else if (loc == "NME_NAXE") {
    scenario <- "Gene-Level Definitive vs. Single-Patient TR Locus Expansion"
    rationale <- "NAXE coding LoF mutations are Definitive/Strong for autosomal recessive mitochondrial neurometabolic encephalopathy, whereas the biallelic promoter GGGCC expansion has only been reported in a single proband to date, yielding Limited locus evidence in criTRia."
  } else if (grepl("^HFG_HOXA13", loc)) {
    scenario <- "Gene-Level Definitive vs. Specific Polyalanine Tract Locus Resolution"
    rationale <- paste0("HOXA13 gene disruptions are Definitive for Hand-Foot-Genital syndrome, but criTRia curates each individual polyalanine tract separately; this specific tract (", d_id, ") currently has Limited published proband/family replication.")
  } else if (loc == "XLID_SOX3") {
    scenario <- "Gene-Level Strong/Definitive vs. Variable/Limited Locus Evidence & Submitter Discordance"
    rationale <- "SOX3 coding alterations are Strong/Definitive for X-linked panhypopituitarism across several groups (G2P, Labcorp, ClinGen Moderate), but polyalanine repeat tract mutations have limited isolated case reports with conflicting interpretations across submitters (Ambry/PanelApp: Limited; G2P: Definitive; Labcorp: Strong)."
  } else if (loc == "TOF_TBX1") {
    scenario <- "Gene-Level Definitive vs. Rare Single-Proband TR Locus"
    rationale <- "TBX1 is Definitively associated with conotruncal heart defects (Ambry: Definitive, Labcorp: Strong), but the specific polyalanine expansion locus was observed in only one ToF proband, classified as Limited in criTRia."
  } else if (loc == "VACTERLX_ZIC3") {
    scenario <- "Gene-Level Definitive/Strong vs. Limited Polyalanine Locus Evidence"
    rationale <- "ZIC3 gene disruptions are Definitive/Strong for X-linked VACTERL association across external groups (G2P: Definitive, PanelApp: Strong, Ambry: Moderate), whereas evidence specifically for the N-terminal polyalanine tract expansion is Limited."
  } else if (loc == "FRA12A_DIP2B") {
    scenario <- "TR Locus Classification Discordance Among Submitters"
    rationale <- "DIP2B CGG repeat expansion is evaluated as Moderate in criTRia with functional and multi-cohort evidence, but external submitters have discordant assessments (Ambry/G2P: Limited; Labcorp: No Known Disease Relationship)."
  } else if (loc == "XDP_TAF1") {
    scenario <- "TR Locus vs. Submitter Outlier Assessment"
    rationale <- "TAF1 SVA retrotransposon hexamer expansion is the Definitive causal mechanism for X-linked dystonia-parkinsonism (criTRia: Definitive, Ambry/PanelApp: Strong, G2P: Definitive), but Labcorp classified it as Limited in an earlier assessment."
  } else if (loc %in% c("FAME1_SAMD12", "FAME2_STARD7")) {
    scenario <- "Emerging TR Mechanism vs. Lagging Clinical Panel Inclusion"
    rationale <- paste0(loc, " noncoding repeat expansions are well-established causes of familial adult myoclonic epilepsy (criTRia: Definitive/Strong, Labcorp: Strong), but received Limited classification in PanelApp due to recent discovery or jurisdictional panel criteria.")
  } else if (loc == "OPDM2_GIPC1") {
    scenario <- "Recently Identified TR Locus vs. PanelApp Limited"
    rationale <- "GIPC1 5' UTR CGG expansion is Definitively established for OPDM2 with strong cohort and functional rescue data (criTRia: Definitive, Labcorp: Strong), whereas PanelApp assigned Limited."
  } else if (loc == "OPDM4_RILPL1") {
    scenario <- "Recently Identified TR Locus vs. External Uncurated/Limited"
    rationale <- "RILPL1 5' UTR CGG expansion is Moderate in criTRia with multiple families and characteristic pathology, but Ambry evaluated it as Limited."
  } else if (loc == "SCA31_BEAN1") {
    scenario <- "Complex Pentanucleotide Insertion vs. Ambry Limited"
    rationale <- "BEAN1 intronic complex pentanucleotide insertion is Definitive for SCA31 (criTRia: Definitive, Labcorp: Strong), whereas Ambry classified it as Limited."
  } else if (loc == "SCA37_DAB1") {
    scenario <- "Noncoding Pentanucleotide Insertion vs. Labcorp Limited"
    rationale <- "DAB1 intronic ATTTC insertion is Strong in criTRia and PanelApp with animal models and segregation, but Labcorp classified it as Limited."
  } else {
    if (critria_is_pos && gencc_has_neg) {
      scenario <- "criTRia Positive vs. External Submitter Negative/Limited"
      rationale <- "criTRia locus curation incorporates recent repeat-specific genetic and experimental evidence, whereas one or more external submitters assigned Limited or No Known."
    } else if (critria_is_neg && gencc_has_pos) {
      scenario <- "criTRia Negative/Limited vs. External Submitter Positive"
      rationale <- "External submitters evaluated the overall gene as Definitive/Strong, whereas the specific TR locus is Limited or Disputed in criTRia."
    } else {
      scenario <- "Discordance Across External Submitters"
      rationale <- "External submitters hold discordant classifications ranging between Definitive/Strong/Moderate and Limited/Contradictory/No Known."
    }
  }
  
  s3b_list[[length(s3b_list) + 1]] <- tibble(
    Gene = g,
    Locus_ID = loc,
    Disease_Phenotype = d_id,
    GenCC_Scores = gencc_scores_str,
    criTRia_Classification = cr_class,
    Total_Score = tot_sc,
    Genetic_Score = gen_sc,
    Experimental_Score = exp_sc,
    Discordance_Scenario = scenario,
    Rationale = rationale
  )
}

s3b_table <- bind_rows(s3b_list)

# Export Supplemental Table S3b
write_tsv(s3b_table, "Supplemental_Table_S3b_matched_discordance.tsv")

cat("\nSuccessfully generated:\n")
cat(" - Supplemental_Table_S3b_matched_discordance.tsv (", nrow(s3b_table), " rows)\n", sep="")
print(s3b_table %>% select(Gene, Locus_ID, Disease_Phenotype, criTRia_Classification, Discordance_Scenario))
