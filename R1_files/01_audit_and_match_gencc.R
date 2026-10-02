# analysis/01_audit_and_match_gencc.R
# Re-evaluates GenCC submissions against criTRia loci by validating BOTH gene_symbol AND disease entity.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

# Create output dir if needed
dir.create("analysis/output", recursive = TRUE, showWarnings = FALSE)

# 1. Load data
critria_file <- list.files(pattern = "Supplemental_File_4_criTRia-curations.*\\.tsv")[1]
if (is.na(critria_file)) stop("Supplemental_File_4 criTRia file not found!")

critria <- read_tsv(critria_file, show_col_types = FALSE)
gencc <- read_tsv("gencc-submissions.tsv", show_col_types = FALSE)

cat("Loaded criTRia loci:", nrow(critria), "\n")
cat("Loaded GenCC submissions:", nrow(gencc), "\n\n")

# Mapping of GenCC submitter titles to standardized groups
submitter_map <- c(
  "Ambry Genetics" = "Ambry",
  "ClinGen" = "ClinGen",
  "G2P" = "G2P",
  "Genomics England PanelApp" = "PanelApp",
  "PanelApp Australia" = "PanelApp",
  "Illumina" = "Illumina",
  "Labcorp Genetics (formerly Invitae)" = "Labcorp",
  "Laboratory for Molecular Medicine" = "Lab MM",
  "Myriad Women's Health" = "Myriad",
  "Orphanet" = "Orphanet"
)

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

# Crosswalk matching logic per Section 3
evaluate_gencc_for_locus <- function(locus_id, gene_sym, disease_id, gencc_row) {
  g_search <- if (gene_sym == "NIPA") "NIPA1" else gene_sym
  if (gencc_row$gene_symbol != g_search) {
    return(list(matched = FALSE, reason = "Gene symbol mismatch"))
  }
  
  d_title <- tolower(gencc_row$disease_title)
  d_orig <- tolower(paste(gencc_row$disease_original_title, gencc_row$submitted_as_disease_name))
  all_text <- paste(d_title, d_orig)
  curie <- gencc_row$disease_curie
  
  # ABCD3 (OPDM5_ABCD3)
  if (locus_id == "OPDM5_ABCD3") {
    if (grepl("oculopharyngodistal|opdm", all_text)) {
      return(list(matched = TRUE, reason = "Matched OPDM5"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-OPDM5 (e.g. peroxisomal/bile acid defect)"))
    }
  }
  
  # AFF2 (FRAXE_AFF2)
  if (locus_id == "FRAXE_AFF2") {
    if (grepl("fraxe|fragile x mental retardation 2|fraxe-associated", all_text)) {
      return(list(matched = TRUE, reason = "Matched FRAXE"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-FRAXE"))
    }
  }
  
  # AFF3 (FRA2A_AFF3)
  if (locus_id == "FRA2A_AFF3") {
    if (grepl("fra2a|kinsship|aff3-related|mesomelic|intellectual disability", all_text)) {
      return(list(matched = TRUE, reason = "Matched FRA2A / AFF3-related"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-AFF3-related phenotype"))
    }
  }
  
  # AR (SBMA_AR)
  if (locus_id == "SBMA_AR") {
    if (grepl("androgen insensitivity|prostate cancer", all_text)) {
      return(list(matched = FALSE, reason = "Excluded non-SBMA (e.g. androgen insensitivity/prostate cancer)"))
    }
    if (grepl("spinal and bulbar muscular atrophy|sbma|kennedy", all_text)) {
      return(list(matched = TRUE, reason = "Matched SBMA"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-SBMA AR phenotype"))
    }
  }
  
  # ARX (EIEE1_ARX & PRTS_ARX)
  if (locus_id == "EIEE1_ARX") {
    if (grepl("lissencephaly|hydranencephaly|proud|partington", all_text)) {
      return(list(matched = FALSE, reason = "Excluded non-EIEE1 ARX phenotype"))
    }
    if (grepl("epilep|infantile spasm|west syndrome|eiee|encephalopathy", all_text)) {
      return(list(matched = TRUE, reason = "Matched EIEE1"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-EIEE1 phenotype"))
    }
  }
  
  if (locus_id == "PRTS_ARX") {
    if (grepl("partington", all_text)) {
      return(list(matched = TRUE, reason = "Matched Partington syndrome"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-Partington ARX phenotype"))
    }
  }
  
  # ATN1 (DRPLA_ATN1)
  if (locus_id == "DRPLA_ATN1") {
    if (grepl("dentatorubral|drpla", all_text)) {
      return(list(matched = TRUE, reason = "Matched DRPLA"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-DRPLA (e.g. CHEDDA / congenital hypotonia)"))
    }
  }
  
  # ATXN1 (SCA1_ATXN1)
  if (locus_id == "SCA1_ATXN1") {
    if (grepl("spinocerebellar ataxia.*1\\b|sca1\\b|spinocerebellar ataxia type 1", all_text)) {
      return(list(matched = TRUE, reason = "Matched SCA1"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-SCA1"))
    }
  }
  
  # ATXN2 (SCA2_ATXN2)
  if (locus_id == "SCA2_ATXN2") {
    if (grepl("spinocerebellar ataxia.*2\\b|sca2\\b|spinocerebellar ataxia type 2", all_text)) {
      return(list(matched = TRUE, reason = "Matched SCA2"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-SCA2"))
    }
  }
  
  # ATXN3 (SCA3_ATXN3)
  if (locus_id == "SCA3_ATXN3") {
    if (grepl("spinocerebellar ataxia.*3|sca3|machado-joseph|mjd", all_text)) {
      return(list(matched = TRUE, reason = "Matched SCA3 / MJD"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-SCA3"))
    }
  }
  
  # ATXN7 (SCA7_ATXN7)
  if (locus_id == "SCA7_ATXN7") {
    if (grepl("spinocerebellar ataxia.*7\\b|sca7\\b|spinocerebellar ataxia type 7", all_text)) {
      return(list(matched = TRUE, reason = "Matched SCA7"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-SCA7"))
    }
  }
  
  # ATXN8OS (SCA8_ATXN8OS)
  if (locus_id == "SCA8_ATXN8OS") {
    if (grepl("spinocerebellar ataxia.*8|sca8", all_text)) {
      return(list(matched = TRUE, reason = "Matched SCA8"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-SCA8"))
    }
  }
  
  # ATXN10 (SCA10_ATXN10)
  if (locus_id == "SCA10_ATXN10") {
    if (grepl("spinocerebellar ataxia.*10|sca10", all_text)) {
      return(list(matched = TRUE, reason = "Matched SCA10"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-SCA10"))
    }
  }
  
  # BEAN1 (SCA31_BEAN1)
  if (locus_id == "SCA31_BEAN1") {
    if (grepl("spinocerebellar ataxia.*31|sca31", all_text)) {
      return(list(matched = TRUE, reason = "Matched SCA31"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-SCA31"))
    }
  }
  
  # C9orf72 (FTDALS1_C9orf72)
  if (locus_id == "FTDALS1_C9orf72") {
    if (grepl("amyotrophic lateral sclerosis|frontotemporal dementia|als-ftd|ftdals", all_text)) {
      return(list(matched = TRUE, reason = "Matched FTDALS1"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-FTDALS1"))
    }
  }
  
  # CACNA1A (SCA6_CACNA1A)
  if (locus_id == "SCA6_CACNA1A") {
    if (grepl("episodic ataxia|hemiplegic migraine|epileptic encephalopathy|torticollis|lennox", all_text)) {
      return(list(matched = FALSE, reason = "Excluded non-SCA6 CACNA1A phenotype (EA2/FHM1/DEE42)"))
    }
    if (grepl("spinocerebellar ataxia.*6|sca6", all_text)) {
      return(list(matched = TRUE, reason = "Matched SCA6"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-SCA6"))
    }
  }
  
  # CBL (JBS_CBL)
  if (locus_id == "JBS_CBL") {
    if (grepl("jacobsen|fra11b", all_text)) {
      return(list(matched = TRUE, reason = "Matched Jacobsen syndrome / FRA11B"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-Jacobsen CBL phenotype (e.g. Noonan-like / JMML / CBL-related)"))
    }
  }
  
  # CNBP (DM2_CNBP)
  if (locus_id == "DM2_CNBP") {
    if (grepl("myotonic dystrophy.*2|dm2|promm", all_text)) {
      return(list(matched = TRUE, reason = "Matched DM2"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-DM2"))
    }
  }
  
  # COMP (EDM1-PSACH_COMP)
  if (locus_id == "EDM1-PSACH_COMP") {
    if (grepl("pseudoachondroplasia|multiple epiphyseal dysplasia|edm1|psach", all_text)) {
      return(list(matched = TRUE, reason = "Matched EDM1/PSACH"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-EDM1/PSACH"))
    }
  }
  
  # CSNK1E (EPM_CSNK1E)
  if (locus_id == "EPM_CSNK1E") {
    if (grepl("progressive myoclon|epm", all_text)) {
      return(list(matched = TRUE, reason = "Matched Progressive myoclonus epilepsy"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-EPM CSNK1E phenotype"))
    }
  }
  
  # CSTB (EPM1_CSTB)
  if (locus_id == "EPM1_CSTB") {
    if (grepl("progressive myoclonus epilepsy 1|epm1|unverricht-lundborg", all_text)) {
      return(list(matched = TRUE, reason = "Matched EPM1"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-EPM1"))
    }
  }
  
  # DAB1 (SCA37_DAB1)
  if (locus_id == "SCA37_DAB1") {
    if (grepl("spinocerebellar ataxia.*37|sca37", all_text)) {
      return(list(matched = TRUE, reason = "Matched SCA37"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-SCA37"))
    }
  }
  
  # DIP2B (FRA12A_DIP2B)
  if (locus_id == "FRA12A_DIP2B") {
    if (grepl("fra12a|intellectual disability.*fra12a|fragile site", all_text)) {
      return(list(matched = TRUE, reason = "Matched FRA12A"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-FRA12A"))
    }
  }
  
  # DMD (DMD_DMD)
  if (locus_id == "DMD_DMD") {
    if (grepl("duchenne|becker|progressive muscular dystrophy", all_text)) {
      return(list(matched = TRUE, reason = "Matched DMD/BMD"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-DMD/BMD"))
    }
  }
  
  # DMPK (DM1_DMPK)
  if (locus_id == "DM1_DMPK") {
    if (grepl("myotonic dystrophy.*1|dm1|steinert", all_text)) {
      return(list(matched = TRUE, reason = "Matched DM1"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-DM1"))
    }
  }
  
  # EIF4A3 (RCPS_EIF4A3)
  if (locus_id == "RCPS_EIF4A3") {
    if (grepl("richieri|rcps", all_text)) {
      return(list(matched = TRUE, reason = "Matched RCPS"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-RCPS"))
    }
  }
  
  # FGF14 (SCA27B_FGF14)
  if (locus_id == "SCA27B_FGF14") {
    if (grepl("27b|gaa-fgf14", all_text)) {
      return(list(matched = TRUE, reason = "Matched SCA27B"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-SCA27B (e.g. SCA27A / SCA27 point mutations)"))
    }
  }
  
  # FMR1 (FXS_FMR1 & FXTAS,POF1_FMR1)
  if (locus_id == "FXS_FMR1") {
    if (grepl("fragile x syndrome|fxs", all_text)) {
      return(list(matched = TRUE, reason = "Matched FXS"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-FXS"))
    }
  }
  
  if (locus_id %in% c("FXTAS,POF1_FMR1", "FXTAS_FMR1", "POF1_FMR1")) {
    if (grepl("tremor/ataxia|fxtas|premature ovarian failure|pof1", all_text)) {
      return(list(matched = TRUE, reason = "Matched FXTAS/POF1"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-FXTAS/POF1"))
    }
  }
  
  # FOXL2 (BPES_FOXL2)
  if (locus_id == "BPES_FOXL2") {
    if (grepl("blepharophimosis|bpes", all_text)) {
      return(list(matched = TRUE, reason = "Matched BPES"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-BPES (e.g. isolated POF3)"))
    }
  }
  
  # FXN (FRDA_FXN)
  if (locus_id == "FRDA_FXN") {
    if (grepl("friedreich|frda", all_text)) {
      return(list(matched = TRUE, reason = "Matched Friedreich ataxia"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-FRDA"))
    }
  }
  
  # GIPC1 (OPDM2_GIPC1)
  if (locus_id == "OPDM2_GIPC1") {
    if (grepl("oculopharyngodistal|opdm", all_text)) {
      return(list(matched = TRUE, reason = "Matched OPDM2"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-OPDM2"))
    }
  }
  
  # GLS (GDPAG_GLS)
  if (locus_id == "GDPAG_GLS") {
    if (grepl("developmental and epileptic encephalopathy|glutaminase deficiency|gdpag", all_text)) {
      return(list(matched = TRUE, reason = "Matched GDPAG / Glutaminase deficiency"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-GDPAG"))
    }
  }
  
  # HOXA13 (HFG-I, HFG-II, HFG-III)
  if (grepl("^HFG", locus_id) || grepl("HOXA13", locus_id)) {
    if (grepl("hand-foot-genital|hfg", all_text)) {
      return(list(matched = TRUE, reason = "Matched Hand-foot-genital syndrome"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-HFG (e.g. Guttmacher)"))
    }
  }
  
  # HOXD13 (SD5_HOXD13)
  if (locus_id == "SD5_HOXD13") {
    if (grepl("synpolydactyly|spd|syndactyly type 5|syndactyly type v\\b|brachydactyly-syndactyly", all_text)) {
      return(list(matched = TRUE, reason = "Matched Synpolydactyly / SD5"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-SD5 (e.g. isolated brachydactyly type D/E)"))
    }
  }
  
  # HTT (HD_HTT)
  if (locus_id == "HD_HTT") {
    if (grepl("huntington disease|\\bhd\\b", all_text)) {
      return(list(matched = TRUE, reason = "Matched Huntington disease"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-HD"))
    }
  }
  
  # JPH3 (HDL2_JPH3)
  if (locus_id == "HDL2_JPH3") {
    if (grepl("huntington disease-like 2|hdl2", all_text)) {
      return(list(matched = TRUE, reason = "Matched HDL2"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-HDL2"))
    }
  }
  
  # LRP12 (OPDM1_LRP12)
  if (locus_id == "OPDM1_LRP12") {
    if (grepl("oculopharyngodistal|opdm", all_text)) {
      return(list(matched = TRUE, reason = "Matched OPDM1"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-OPDM1"))
    }
  }
  
  # MARCHF6 (FAME3_MARCHF6)
  if (locus_id == "FAME3_MARCHF6") {
    if (grepl("myoclonic epilepsy|fame|bafme", all_text)) {
      return(list(matched = TRUE, reason = "Matched FAME3"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-FAME3"))
    }
  }
  
  # MIR7-2 (CHNG3_MIR7-2)
  if (locus_id == "CHNG3_MIR7-2") {
    if (grepl("congenital heart disease, non-syndromic, 3|chng3|thyroid|hypothyroidism", all_text)) {
      return(list(matched = TRUE, reason = "Matched CHNG3"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-CHNG3"))
    }
  }
  
  # MUC1 (ADTKD_MUC1)
  if (locus_id == "ADTKD_MUC1") {
    if (grepl("tubulointerstitial kidney disease|adtkd|medullary cystic kidney", all_text)) {
      return(list(matched = TRUE, reason = "Matched ADTKD"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-ADTKD"))
    }
  }
  
  # NAXE (NME_NAXE)
  if (locus_id == "NME_NAXE") {
    if (grepl("encephalopathy|pebat|nad\\(p\\)hx|mitochondrial", all_text)) {
      return(list(matched = TRUE, reason = "Matched NME / mitochondrial encephalopathy"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-NME"))
    }
  }
  
  # NIPA1 (ALS1_NIPA1)
  if (locus_id %in% c("ALS1_NIPA1", "ALS1_NIPA")) {
    if (grepl("spastic paraplegia|spg6", all_text)) {
      return(list(matched = FALSE, reason = "Excluded SPG6 (non-ALS phenotype)"))
    }
    if (grepl("amyotrophic lateral sclerosis|als", all_text)) {
      return(list(matched = TRUE, reason = "Matched ALS1"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-ALS"))
    }
  }
  
  # NOP56 (SCA36_NOP56)
  if (locus_id == "SCA36_NOP56") {
    if (grepl("spinocerebellar ataxia.*36|sca36|asidan", all_text)) {
      return(list(matched = TRUE, reason = "Matched SCA36"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-SCA36"))
    }
  }
  
  # NOTCH2NLC (NIID_NOTCH2NLC)
  if (locus_id == "NIID_NOTCH2NLC") {
    if (grepl("neuronal intranuclear inclusion|niid", all_text)) {
      return(list(matched = TRUE, reason = "Matched NIID"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-NIID"))
    }
  }
  
  # NUTM2B-AS1 (OPML1_NUTM2B_AS1)
  if (grepl("OPML1", locus_id) || grepl("NUTM2B", locus_id)) {
    if (grepl("oculopharyngeal muscular dystrophy-like|opml", all_text)) {
      return(list(matched = TRUE, reason = "Matched OPML1"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-OPML1"))
    }
  }
  
  # PABPN1 (OPMD_PABPN1)
  if (locus_id == "OPMD_PABPN1") {
    if (grepl("oculopharyngeal muscular dystrophy|opmd", all_text)) {
      return(list(matched = TRUE, reason = "Matched OPMD"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-OPMD"))
    }
  }
  
  # PLIN4 (MRUPAV_PLIN4)
  if (locus_id == "MRUPAV_PLIN4") {
    if (grepl("rimmed|vacuol|mrupav|neuromyopathy|myopathy", all_text)) {
      return(list(matched = TRUE, reason = "Matched MRUPAV"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-MRUPAV"))
    }
  }
  
  # POLG (CPEO_POLG)
  if (locus_id == "CPEO_POLG") {
    if (grepl("progressive external ophthalmoplegia|cpeo|peob", all_text)) {
      return(list(matched = TRUE, reason = "Matched CPEO"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-CPEO POLG phenotype (e.g. Leigh / SANDO / MNGIE / depletion syndrome 4a)"))
    }
  }
  
  # PPP2R2B (SCA12_PPP2R2B)
  if (locus_id == "SCA12_PPP2R2B") {
    if (grepl("spinocerebellar ataxia.*12|sca12", all_text)) {
      return(list(matched = TRUE, reason = "Matched SCA12"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-SCA12"))
    }
  }
  
  # PRDM12 (HSAN-VIII_PRDM12)
  if (locus_id == "HSAN-VIII_PRDM12") {
    if (grepl("hereditary sensory|hsan|insensitivity to pain", all_text)) {
      return(list(matched = TRUE, reason = "Matched HSAN-VIII / pain insensitivity"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-HSAN-VIII"))
    }
  }
  
  # PRNP (CJD_PRNP)
  if (locus_id == "CJD_PRNP") {
    if (grepl("creutzfeldt|cjd|gerstmann|gss|prion|huntington disease-like 1|insomnia", all_text)) {
      return(list(matched = TRUE, reason = "Matched Prion disease / CJD / GSS / HDL1"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-CJD/prion phenotype"))
    }
  }
  
  # RAI1 (FAME8_RAI1)
  if (locus_id == "FAME8_RAI1") {
    if (grepl("myoclonic epilepsy|fame8", all_text)) {
      return(list(matched = TRUE, reason = "Matched FAME8"))
    } else {
      return(list(matched = FALSE, reason = "Excluded Smith-Magenis / Potocki-Lupski (non-FAME8)"))
    }
  }
  
  # RAPGEF2 (FAME7_RAPGEF2)
  if (locus_id == "FAME7_RAPGEF2") {
    if (grepl("myoclonic epilepsy|fame7|cortical myoclonic", all_text)) {
      return(list(matched = TRUE, reason = "Matched FAME7"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-FAME7 (e.g. ALS)"))
    }
  }
  
  # RFC1 (CANVAS_RFC1)
  if (locus_id == "CANVAS_RFC1") {
    if (grepl("canvas|cerebellar ataxia, neuropathy, vestibular", all_text)) {
      return(list(matched = TRUE, reason = "Matched CANVAS"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-CANVAS"))
    }
  }
  
  # RILPL1 (OPDM4_RILPL1)
  if (locus_id == "OPDM4_RILPL1") {
    if (grepl("oculopharyngodistal|opdm4", all_text)) {
      return(list(matched = TRUE, reason = "Matched OPDM4"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-OPDM4"))
    }
  }
  
  # RUNX2 (CCD_RUNX2)
  if (locus_id == "CCD_RUNX2") {
    if (grepl("cleidocranial|ccd", all_text)) {
      return(list(matched = TRUE, reason = "Matched Cleidocranial dysplasia"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-CCD (e.g. metaphyseal dysplasia)"))
    }
  }
  
  # SAMD12 (FAME1_SAMD12)
  if (locus_id == "FAME1_SAMD12") {
    if (grepl("myoclonic epilepsy|fame1|bafme1", all_text)) {
      return(list(matched = TRUE, reason = "Matched FAME1"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-FAME1"))
    }
  }
  
  # SOX3 (XLID_SOX3)
  if (locus_id %in% c("XLID_SOX3", "XLID, PHPX_SOX3")) {
    if (grepl("sex reversal|hypertrichosis|septooptic", all_text)) {
      return(list(matched = FALSE, reason = "Excluded SRXX3 / hypertrichosis / septooptic dysplasia"))
    }
    if (grepl("panhypopituitarism|pituitary hormone deficiency|growth hormone deficiency", all_text)) {
      return(list(matched = TRUE, reason = "Matched Panhypopituitarism / PHPX"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-PHPX SOX3 phenotype"))
    }
  }
  
  # STARD7 (FAME2_STARD7)
  if (locus_id == "FAME2_STARD7") {
    if (grepl("myoclonic epilepsy|fame2", all_text)) {
      return(list(matched = TRUE, reason = "Matched FAME2"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-FAME2"))
    }
  }
  
  # TAF1 (XDP_TAF1)
  if (locus_id == "XDP_TAF1") {
    if (grepl("dystonia-parkinsonism|xdp|lubag", all_text)) {
      return(list(matched = TRUE, reason = "Matched XDP"))
    } else {
      return(list(matched = FALSE, reason = "Excluded MRX33 / intellectual disability (non-XDP)"))
    }
  }
  
  # TBP (SCA17_TBP)
  if (locus_id == "SCA17_TBP") {
    if (grepl("spinocerebellar ataxia.*17|sca17", all_text)) {
      return(list(matched = TRUE, reason = "Matched SCA17"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-SCA17"))
    }
  }
  
  # TBX1 (TOF_TBX1)
  if (locus_id == "TOF_TBX1") {
    if (grepl("22q11|digeorge|velocardiofacial", all_text)) {
      return(list(matched = FALSE, reason = "Excluded 22q11.2 deletion / DiGeorge / VCFS"))
    }
    if (grepl("tetralogy of fallot|conotruncal heart malformations|tof", all_text)) {
      return(list(matched = TRUE, reason = "Matched Tetralogy of Fallot / Conotruncal heart malformations"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-TOF TBX1 phenotype"))
    }
  }
  
  # TCF4 (FECD3_TCF4)
  if (locus_id == "FECD3_TCF4") {
    if (grepl("fuchs|corneal dystrophy|fecd", all_text)) {
      return(list(matched = TRUE, reason = "Matched Fuchs endothelial corneal dystrophy 3"))
    } else {
      return(list(matched = FALSE, reason = "Excluded Pitt-Hopkins / non-FECD"))
    }
  }
  
  # THAP11 (SCA51_THAP11)
  if (locus_id == "SCA51_THAP11") {
    if (grepl("spinocerebellar ataxia.*51|sca51", all_text)) {
      return(list(matched = TRUE, reason = "Matched SCA51"))
    } else {
      return(list(matched = FALSE, reason = "Excluded Cobalamin metabolism defect / MMA"))
    }
  }
  
  # TNRC6A (FAME6_TNRC6A)
  if (locus_id == "FAME6_TNRC6A") {
    if (grepl("myoclonic epilepsy|fame6", all_text)) {
      return(list(matched = TRUE, reason = "Matched FAME6"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-FAME6"))
    }
  }
  
  # TYMS (CPUM_TYMS)
  if (locus_id == "CPUM_TYMS") {
    if (grepl("corneal dystrophy|cpum|melanosis", all_text)) {
      return(list(matched = TRUE, reason = "Matched CPUM"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-CPUM"))
    }
  }
  
  # VWA1 (HMNR7_VWA1)
  if (locus_id == "HMNR7_VWA1") {
    if (grepl("hereditary motor neuropathy|hmnr7|neuronopathy, distal hereditary motor", all_text)) {
      return(list(matched = TRUE, reason = "Matched HMNR7"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-HMNR7"))
    }
  }
  
  # XYLT1 (DBQD2_XYLT1)
  if (locus_id == "DBQD2_XYLT1") {
    if (grepl("desbuquois|dbqd", all_text)) {
      return(list(matched = TRUE, reason = "Matched Desbuquois dysplasia 2"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-DBQD2"))
    }
  }
  
  # YEATS2 (FAME4_YEATS2)
  if (locus_id == "FAME4_YEATS2") {
    if (grepl("myoclonic epilepsy|fame4", all_text)) {
      return(list(matched = TRUE, reason = "Matched FAME4"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-FAME4"))
    }
  }
  
  # ZFHX3 (SCA4_ZFHX3)
  if (locus_id == "SCA4_ZFHX3") {
    if (grepl("spinocerebellar ataxia.*4\\b|sca4\\b|spinocerebellar ataxia type 4", all_text)) {
      return(list(matched = TRUE, reason = "Matched SCA4"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-SCA4 (e.g. complex NDD / epilepsy)"))
    }
  }
  
  # ZIC2 (HPE5_ZIC2)
  if (locus_id == "HPE5_ZIC2") {
    if (grepl("holoprosencephaly.*5|hpe5|holoprosencephaly", all_text)) {
      return(list(matched = TRUE, reason = "Matched HPE5"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-HPE5"))
    }
  }
  
  # ZIC3 (VACTERLX_ZIC3)
  if (locus_id == "VACTERLX_ZIC3") {
    if (grepl("heterotaxy", all_text)) {
      return(list(matched = FALSE, reason = "Excluded Heterotaxy 1 (HTX1)"))
    }
    if (grepl("vacterl", all_text)) {
      return(list(matched = TRUE, reason = "Matched VACTERLX"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-VACTERLX"))
    }
  }
  
  # ZNF713 (FRA7A_ZNF713)
  if (locus_id == "FRA7A_ZNF713") {
    if (grepl("fra7a|fragile site 7a", all_text)) {
      return(list(matched = TRUE, reason = "Matched FRA7A"))
    } else {
      return(list(matched = FALSE, reason = "Excluded non-FRA7A (e.g. autism)"))
    }
  }
  
  return(list(matched = FALSE, reason = "No specific match rule satisfied"))
}

# Iterate over all loci
all_matches <- list()
all_exclusions <- list()

for (i in 1:nrow(critria)) {
  loc_id <- critria$Locus_ID[i]
  gene_sym <- critria$Gene[i]
  dis_id <- critria$Disease_ID[i]
  
  g_search <- if (gene_sym == "NIPA") "NIPA1" else gene_sym
  gencc_sub <- gencc %>% filter(gene_symbol == g_search)
  
  if (nrow(gencc_sub) == 0) next
  
  for (j in 1:nrow(gencc_sub)) {
    grow <- gencc_sub[j, ]
    res <- evaluate_gencc_for_locus(loc_id, gene_sym, dis_id, grow)
    sub_title <- grow$submitter_title
    sub_grp <- if (sub_title %in% names(submitter_map)) submitter_map[[sub_title]] else sub_title
    
    if (res$matched) {
      all_matches[[length(all_matches) + 1]] <- tibble(
        Locus_ID = loc_id,
        Gene = gene_sym,
        Disease_ID = dis_id,
        Submitter = sub_grp,
        Classification_Raw = grow$classification_title,
        Score = standardize_class(grow$classification_title),
        Date = grow$submitted_as_date,
        Disease_Title_GenCC = grow$disease_title,
        Disease_Curie_GenCC = grow$disease_curie,
        Submitter_Raw = sub_title
      )
    } else {
      all_exclusions[[length(all_exclusions) + 1]] <- tibble(
        Gene = gene_sym,
        Locus_ID = loc_id,
        Submitter = sub_title,
        Disease_Title_GenCC = grow$disease_title,
        Reason_for_Exclusion = res$reason
      )
    }
  }
}

matches_df <- bind_rows(all_matches)
exclusions_df <- bind_rows(all_exclusions)

# Step 3: For submitters with multiple submissions for the same matched disease entity, take the most recent
# Note: Deduplicate by (Locus_ID, Submitter) taking the latest Date
matches_dedup <- matches_df %>%
  arrange(desc(Date)) %>%
  group_by(Locus_ID, Submitter) %>%
  slice(1) %>%
  ungroup()

# Required columns for matched_gencc_curations.tsv:
# Locus_ID, Gene, Disease_ID, Submitter, Score, Date, Disease_Title_GenCC, Disease_Curie_GenCC
matched_out <- matches_dedup %>%
  select(Locus_ID, Gene, Disease_ID, Submitter, Score, Date, Disease_Title_GenCC, Disease_Curie_GenCC)

# Required columns for excluded_mismatched_curations.tsv:
# Gene, Locus_ID, Submitter, Disease_Title_GenCC, Reason_for_Exclusion
excluded_out <- exclusions_df %>%
  select(Gene, Locus_ID, Submitter, Disease_Title_GenCC, Reason_for_Exclusion)

# Write output TSVs
write_tsv(matched_out, "matched_gencc_curations.tsv")
write_tsv(excluded_out, "excluded_mismatched_curations.tsv")

cat("Successfully generated:\n")
cat(" - matched_gencc_curations.tsv (", nrow(matched_out), " rows)\n", sep="")
cat(" - excluded_mismatched_curations.tsv (", nrow(excluded_out), " rows)\n", sep="")

# Summary statistics
cat("\nSummary Statistics:\n")
cat("Total GenCC gene-level candidate submissions evaluated:", nrow(matches_df) + nrow(exclusions_df), "\n")
cat("Total submissions matched (prior to submitter dedup):", nrow(matches_df), "\n")
cat("Total submissions matched (deduplicated per submitter):", nrow(matched_out), "\n")
cat("Total submissions excluded as disease mismatches:", nrow(excluded_out), "\n")

excluded_loci <- unique(excluded_out$Locus_ID)
cat("\nLoci with excluded GenCC disease mismatches (", length(excluded_loci), " loci):\n", sep="")
cat(paste(excluded_loci, collapse=", "), "\n")
