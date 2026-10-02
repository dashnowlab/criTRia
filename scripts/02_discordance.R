# 02_discordance.R
# Combines criTRia and matched GenCC classifications for the criTRia-curated loci
# into a wide matrix, flags discordance, and builds the discordant loci table.
# Run from the repository root (see run_all.sh).
# Inputs (from scripts/01_download_and_match.py):
#   paper/supp3_dataset.tsv            criTRia and matched GenCC classifications
#   paper/supp4_criTRia_curations.tsv  criTRia curations, with Locus_ID
# Outputs:
#   data/processed/discordance_matrix.tsv
#   paper/suppX_discordant_loci.tsv    Supplemental File X

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})
source("scripts/discordance.R")

# 1. Load data
critria <- read_tsv("paper/supp4_criTRia_curations.tsv", show_col_types = FALSE)

# criTRia and GenCC classifications for the criTRia-curated loci
dataset <- read_tsv("paper/supp3_dataset.tsv", show_col_types = FALSE) %>%
  filter(Locus_ID %in% critria$Locus_ID)

matched_filtered <- dataset %>%
  filter(Group != "criTRia") %>%
  transmute(Locus_ID, Submitter = Group, Score = categorical_score)

# Submitter column order
target_submitters <- c("criTRia", "Labcorp", "Orphanet", "Ambry", "PanelApp", "G2P", "ClinGen", "Lab MM", "Illumina", "Myriad")

# 2. Pivot to a wide matrix, one row per locus in curation order
matrix_wide <- dataset %>%
  pivot_wider(id_cols = Locus_ID, names_from = Group, values_from = categorical_score) %>%
  left_join(critria %>% select(Locus_ID, Gene, Disease_ID), by = "Locus_ID") %>%
  arrange(match(Locus_ID, critria$Locus_ID))

# Fill uncurated cells
for (col in target_submitters) {
  if (!col %in% colnames(matrix_wide)) {
    matrix_wide[[col]] <- "Uncurated"
  } else {
    matrix_wide[[col]] <- ifelse(is.na(matrix_wide[[col]]), "Uncurated", matrix_wide[[col]])
  }
}

# 3. Add discordance flags (see ../discordance.R)
matrix_wide <- matrix_wide %>%
  select(Locus_ID, Gene, Disease_ID, all_of(target_submitters)) %>%
  left_join(
    discordance_flags(dataset) %>%
      rename(Conflict_vs_criTRia = vs_criTRia, Conflict_within_GenCC = within_GenCC),
    by = "Locus_ID"
  )

write_tsv(matrix_wide, "data/processed/discordance_matrix.tsv")

cat("Discordance matrix computed:\n")
cat(" - Total loci:", nrow(matrix_wide), "\n")
cat(" - Discordant vs criTRia:", sum(matrix_wide$Conflict_vs_criTRia), "\n")
cat(" - Discordant within GenCC:", sum(matrix_wide$Conflict_within_GenCC), "\n")

# 4. Build the discordant loci table (Supplemental File X)

conflict_loci_df <- matrix_wide %>% filter(Conflict_vs_criTRia | Conflict_within_GenCC)

s3b_list <- list()

for (i in 1:nrow(conflict_loci_df)) {
  loc <- conflict_loci_df$Locus_ID[i]
  g <- conflict_loci_df$Gene[i]
  d_id <- conflict_loci_df$Disease_ID[i]
  
  # criTRia info
  cr_info <- critria %>% filter(Locus_ID == loc)
  cr_class <- conflict_loci_df$criTRia[i]
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
    rationale <- "GenCC submitters disagree on SOX3 for X-linked hypopituitarism with or without intellectual disability (G2P: Definitive; Labcorp: Strong; ClinGen: Moderate; Ambry and PanelApp: Limited), whereas criTRia classifies the polyalanine repeat expansion as Limited (5.0 points) based on a small number of case reports."
  } else if (loc == "TOF_TBX1") {
    scenario <- "Gene-Level Definitive vs. Rare Single-Proband TR Locus"
    rationale <- "TBX1 is Definitively associated with conotruncal heart defects (Ambry: Definitive, Labcorp: Strong), but the specific polyalanine expansion locus was observed in only one ToF proband, classified as Limited in criTRia."
  } else if (loc == "VACTERLX_ZIC3") {
    scenario <- "Gene-Level Definitive/Strong vs. Limited Polyalanine Locus Evidence"
    rationale <- "ZIC3 gene disruptions are Definitive/Strong for X-linked VACTERL association across external groups (G2P: Definitive, PanelApp: Strong, Ambry: Moderate), whereas evidence specifically for the N-terminal polyalanine tract expansion is Limited."
  } else if (loc == "XDP_TAF1") {
    scenario <- "TR Locus vs. Submitter Outlier Assessment"
    rationale <- "The TAF1 SVA retrotransposon hexamer expansion is the established cause of X-linked dystonia-parkinsonism (criTRia: Definitive; Ambry and PanelApp: Strong), but Labcorp classified the relationship as Limited."
  } else if (loc == "FAME1_SAMD12") {
    scenario <- "Emerging TR Mechanism vs. Lagging Clinical Panel Inclusion"
    rationale <- "SAMD12 intronic repeat expansions are an established cause of familial adult myoclonic epilepsy type 1 (criTRia: Definitive; Labcorp: Strong), but PanelApp classified the relationship as Limited."
  } else if (loc == "FAME2_STARD7") {
    scenario <- "Emerging TR Mechanism vs. Lagging Clinical Panel Inclusion"
    rationale <- "STARD7 intronic repeat expansions are an established cause of familial adult myoclonic epilepsy type 2 (criTRia: Strong), but PanelApp, the only GenCC submitter, classified the relationship as Limited."
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
    Conflict_vs_criTRia = conflict_loci_df$Conflict_vs_criTRia[i],
    Conflict_within_GenCC = conflict_loci_df$Conflict_within_GenCC[i],
    Total_Score = tot_sc,
    Genetic_Score = gen_sc,
    Experimental_Score = exp_sc,
    Discordance_Scenario = scenario,
    Rationale = rationale
  )
}

s3b_table <- bind_rows(s3b_list)

# Export the discordant loci table
write_tsv(s3b_table, "paper/suppX_discordant_loci.tsv")

cat("\nSuccessfully generated:\n")
cat(" - paper/suppX_discordant_loci.tsv (", nrow(s3b_table), " rows)\n", sep="")
print(s3b_table %>% select(Gene, Locus_ID, Disease_Phenotype, criTRia_Classification, Discordance_Scenario))
