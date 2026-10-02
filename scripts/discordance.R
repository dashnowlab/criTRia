# Discordance flags, shared by scripts/02_discordance.R and scripts/03_figures.R.
#
# A high score is Definitive, Strong or Moderate; a low score is Limited,
# Contradictory or No Known. Supportive (Orphanet) counts as neither.
#   vs_criTRia:   criTRia is on the opposite side from at least one GenCC group.
#   within_GenCC: GenCC groups disagree among themselves.

high_evidence <- c("Definitive", "Strong", "Moderate")
low_evidence  <- c("Limited", "Contradictory", "No Known")

# data: one row per locus and group, with columns Locus_ID, Group and
# categorical_score. Returns one row per locus with the two flags.
discordance_flags <- function(data) {
  data %>%
    group_by(Locus_ID) %>%
    summarise(
      critria_high = any(Group == "criTRia" & categorical_score %in% high_evidence),
      critria_low  = any(Group == "criTRia" & categorical_score %in% low_evidence),
      gencc_high   = any(Group != "criTRia" & categorical_score %in% high_evidence),
      gencc_low    = any(Group != "criTRia" & categorical_score %in% low_evidence),
      .groups = "drop"
    ) %>%
    transmute(
      Locus_ID,
      vs_criTRia   = (critria_high & gencc_low) | (critria_low & gencc_high),
      within_GenCC = gencc_high & gencc_low
    )
}
