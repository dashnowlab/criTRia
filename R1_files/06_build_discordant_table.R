# 06_build_discordant_table.R
# Builds Supplemental_File_6_discordant_loci_R2.tsv
# from Supplemental_Table_S3b_matched_discordance.tsv and Supplemental_File_4_criTRia-curations_R1.tsv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

# 1. Load inputs
critria_files <- list.files(pattern = "Supplemental_File_4_criTRia-curations.*\\.tsv")
if (length(critria_files) == 0) stop("Supplemental_File_4 criTRia file not found in current directory!")
critria_file <- critria_files[1]

s3b_file <- "Supplemental_Table_S3b_matched_discordance.tsv"
if (!file.exists(s3b_file)) stop("Supplemental_Table_S3b_matched_discordance.tsv not found in current directory!")

critria <- read_tsv(critria_file, show_col_types = FALSE)
s3b <- read_tsv(s3b_file, show_col_types = FALSE)

# 2. Exact 17 Loci Order
locus_order <- c(
  "DMD_DMD", "CPEO_POLG",
  "HFG_HOXA13-I", "HFG_HOXA13-II", "HFG_HOXA13-III",
  "NME_NAXE", "TOF_TBX1", "VACTERLX_ZIC3",
  "SCA31_BEAN1", "FAME1_SAMD12", "FAME2_STARD7", "OPDM2_GIPC1", "OPDM4_RILPL1", "SCA37_DAB1",
  "FRA12A_DIP2B", "XLID_SOX3", "XDP_TAF1"
)

category_map <- c(
  "DMD_DMD" = "Scenario 1: Locus vs. Gene Discordance",
  "CPEO_POLG" = "Scenario 1: Locus vs. Gene Discordance",
  "HFG_HOXA13-I" = "Scenario 2: Intragenic Tract Resolution",
  "HFG_HOXA13-II" = "Scenario 2: Intragenic Tract Resolution",
  "HFG_HOXA13-III" = "Scenario 2: Intragenic Tract Resolution",
  "NME_NAXE" = "Scenario 3: Emerging Single-Proband / Low-Penetrance TR Expansion",
  "TOF_TBX1" = "Scenario 3: Emerging Single-Proband / Low-Penetrance TR Expansion",
  "VACTERLX_ZIC3" = "Scenario 3: Emerging Single-Proband / Low-Penetrance TR Expansion",
  "SCA31_BEAN1" = "Scenario 4: Diagnostic Panel Curation Lag",
  "FAME1_SAMD12" = "Scenario 4: Diagnostic Panel Curation Lag",
  "FAME2_STARD7" = "Scenario 4: Diagnostic Panel Curation Lag",
  "OPDM2_GIPC1" = "Scenario 4: Diagnostic Panel Curation Lag",
  "OPDM4_RILPL1" = "Scenario 4: Diagnostic Panel Curation Lag",
  "SCA37_DAB1" = "Scenario 4: Diagnostic Panel Curation Lag",
  "FRA12A_DIP2B" = "Scenario 5: Inter-Submitter Discordance",
  "XLID_SOX3" = "Scenario 5: Inter-Submitter Discordance",
  "XDP_TAF1" = "Scenario 5: Inter-Submitter Discordance"
)

detailed_descriptions <- c(
  "DMD_DMD" = "ClinGen curates DMD as Definitive for Duchenne muscular dystrophy based on extensive loss-of-function variants across the coding region. In contrast, criTRia specifically curates the 5' UTR polyalanine repeat tract expansion, which lacks robust independent segregating evidence (criTRia score: Limited, 2.0 pts). This illustrates how gene-level curation masks variant-class specific validity.",
  "CPEO_POLG" = "PanelApp and ClinGen curate POLG as Definitive for progressive external ophthalmoplegia and mitochondrial DNA depletion syndromes based on extensive missense and truncating SNVs. criTRia evaluates the exon 2 CAG repeat expansion specifically, which shows contradictory/disputed evidence in recent population cohorts (criTRia score: Limited/Disputed).",
  "HFG_HOXA13-I" = "HOXA13 contains three distinct polyalanine tracts associated with Hand-Foot-Genital Syndrome (HFGS). criTRia curates each tract independently, rating Tract I as Definitive (16.0 pts, fully penetrant expansions). ClinGen and PanelApp evaluate HOXA13 only as a whole gene entity (Definitive), obscuring tract-specific pathogenic thresholds.",
  "HFG_HOXA13-II" = "criTRia evaluates HOXA13 Tract II polyalanine expansions as Moderate (10.0 pts), reflecting intermediate clinical evidence and variable penetrance. Whole-gene ClinGen curation classifies HOXA13 as Definitive, failing to distinguish the intermediate evidence of Tract II from Tract I.",
  "HFG_HOXA13-III" = "criTRia evaluates HOXA13 Tract III expansions as Limited (1.0 pt) due to scarce proband reports and borderline allele length shifts. Gene-level Definitive classification by external submitters overlooks the insufficient evidence for Tract III expansions.",
  "NME_NAXE" = "PanelApp curates NAXE as Definitive/Green for neurometabolic encephalopathy based on biallelic LoF/missense mutations. criTRia evaluates an isolated intragenic polyalanine repeat expansion reported in a single family (criTRia score: Limited, 1.0 pt), distinguishing the emerging TR mechanism from validated null alleles.",
  "TOF_TBX1" = "ClinGen and PanelApp curate TBX1 as Definitive for conotruncal heart defects and 22q11.2 deletion syndrome features. criTRia curates the specific intragenic polyalanine tract expansion polymorphism, which is evaluated as Limited (1.0 pt) due to incomplete penetrance and variable expressivity in isolated TOF.",
  "VACTERLX_ZIC3" = "PanelApp classifies ZIC3 as Definitive for X-linked heterotaxy and congenital heart defects. criTRia scores the specific polyalanine expansion as Limited (2.0 pts) because polyalanine expansions represent rare, candidate alleles relative to established loss-of-function variants.",
  "SCA31_BEAN1" = "criTRia curates the BEAN1 intronic pentanucleotide (TGGAA)n insertion/expansion as Definitive (16.0 pts) for Spinocerebellar Ataxia 31 based on international familial segregations. In contrast, historical gene-level diagnostic panels (PanelApp Australia / G2P) lagged behind, classifying BEAN1 as No Known / Uncurated.",
  "FAME1_SAMD12" = "criTRia classifies SAMD12 intronic TTTCA/TTTTA repeat expansions as Definitive (16.0 pts) for Familial Adult Myoclonic Epilepsy 1 (FAME1). Several routine commercial panels historically lacked non-coding repeat expansion testing, leaving SAMD12 uncurated or disputed on gene panels.",
  "FAME2_STARD7" = "criTRia classifies STARD7 intronic ATTTC repeat expansions as Definitive (16.0 pts) for FAME2. External gene-centric submitters (e.g., G2P/PanelApp) show curation lag, omitting non-coding repeat mechanisms from standard exome panels.",
  "OPDM2_GIPC1" = "criTRia curates the 5' UTR GGC repeat expansion in GIPC1 as Definitive (16.0 pts) for Oculopharyngodistal Myopathy 2 (OPDM2). ClinGen and older diagnostic panels have not established a gene-disease curation for GIPC1, resulting in uncurated/no known status in GenCC.",
  "OPDM4_RILPL1" = "criTRia curates the 5' UTR CGG/CCG repeat expansion in RILPL1 as Definitive (15.0 pts) for OPDM4. External diagnostic submitters lack specialized repeat-expansion curations for RILPL1, representing diagnostic panel curation lag.",
  "SCA37_DAB1" = "criTRia curates the DAB1 5' UTR pentanucleotide (ATTTC)n insertion as Definitive (16.0 pts) for Spinocerebellar Ataxia 37. External panels (e.g. PanelApp/G2P) have historically classified DAB1 as Limited or No Known on gene-level ataxia panels.",
  "FRA12A_DIP2B" = "Genomics England PanelApp classifies DIP2B as Definitive for fragile site 12A / intellectual disability, whereas G2P classifies DIP2B as Limited / Disputed. criTRia independently curates the 5' UTR CGG expansion as Definitive (16.0 pts).",
  "XLID_SOX3" = "Genomics England PanelApp curates SOX3 as Definitive for intellectual disability, whereas ClinGen curates SOX3 as Limited for isolated growth hormone deficiency / hypopituitarism. criTRia scores the polyalanine expansion specifically as Definitive (16.0 pts) for syndromic panhypopituitarism.",
  "XDP_TAF1" = "PanelApp and Orphanet evaluate TAF1 as Definitive for X-linked Dystonia-Parkinsonism (Lubag), whereas other submitters historically flagged TAF1 SVA retrotransposon/hexamer repeat insertions as Limited or Disputed on standard sequencing panels. criTRia classifies XDP_TAF1 as Definitive (16.0 pts)."
)

# Extract and build rows
output_list <- list()

for (loc in locus_order) {
  s3b_row <- s3b %>% filter(Locus_ID == loc)
  cr_row <- critria %>% filter(Locus_ID == loc)
  
  gene_val <- if (nrow(s3b_row) > 0 && !is.na(s3b_row$Gene)) s3b_row$Gene[1] else cr_row$Gene[1]
  dis_name <- if (nrow(s3b_row) > 0 && "Disease_Name" %in% colnames(s3b_row)) s3b_row$Disease_Name[1] else cr_row$Disease_Name[1]
  dis_id <- if (nrow(s3b_row) > 0 && "Disease_ID" %in% colnames(s3b_row)) s3b_row$Disease_ID[1] else cr_row$Disease_ID[1]
  inh_val <- if (nrow(s3b_row) > 0 && "Inheritance" %in% colnames(s3b_row)) s3b_row$Inheritance[1] else cr_row$inheritance[1]
  cr_score <- if (nrow(s3b_row) > 0 && "criTRia_Score" %in% colnames(s3b_row)) s3b_row$criTRia_Score[1] else cr_row$total_score[1]
  cr_class <- if (nrow(s3b_row) > 0 && "criTRia_Classification" %in% colnames(s3b_row)) s3b_row$criTRia_Classification[1] else cr_row$classification[1]
  pos_sub <- if (nrow(s3b_row) > 0 && "Positive_Submitters" %in% colnames(s3b_row)) s3b_row$Positive_Submitters[1] else ""
  neg_sub <- if (nrow(s3b_row) > 0 && "Negative_Submitters" %in% colnames(s3b_row)) s3b_row$Negative_Submitters[1] else ""
  all_sub <- if (nrow(s3b_row) > 0 && "All_Evaluated_Submitters" %in% colnames(s3b_row)) s3b_row$All_Evaluated_Submitters[1] else ""
  
  output_list[[length(output_list) + 1]] <- tibble(
    Locus_ID = loc,
    Gene = gene_val,
    Disease_Name = dis_name,
    Disease_ID = dis_id,
    Inheritance = inh_val,
    criTRia_Score = cr_score,
    criTRia_Classification = cr_class,
    Discordance_Category = category_map[[loc]],
    Positive_Submitters = pos_sub,
    Negative_Submitters = neg_sub,
    All_Evaluated_Submitters = all_sub,
    Scientific_Mechanism_and_Clinical_Relevance = detailed_descriptions[[loc]]
  )
}

output_df <- bind_rows(output_list)

tsv_out <- "Supplemental_File_6_discordant_loci_R2.tsv"
write_tsv(output_df, tsv_out)

cat("Successfully generated:", tsv_out, "(", nrow(output_df), "rows)\n")
