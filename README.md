# criTRia

criTRia (pronounce criteria) is a formalized scoring framework designed to accurately evaluate TR-disease relationships. criTRia builds on current best-practice scoring procedures developed by the Clinical Genome Resource (ClinGen) while (re)curating based on locus-level (rather than gene-level) classification, introducing TR-specific evidence categories and reweighted scoring. This repository contains the data and code for the criTRia manuscript, and always reflects the latest version of the paper. Earlier versions are available from the git history.

Preprint: Weiner MA, Hiatt L, Ajuyah P, Aliyev E, Dashnow H. criTRia: A Classification System and Evidence Criteria for Tandem Repeat Locus-Disease Relationships. medRxiv (2026). https://doi.org/10.64898/2026.07.04.26357279

The current version of the criTRia curations is on the [STRchive criTRia page](https://strchive.org/critria), which lists each locus-disease classification with its score and links to the criTRia SOP. The curations are maintained in the STRchive GitHub repository as [`data/criTRia-curations.tsv`](https://github.com/dashnowlab/STRchive/blob/main/data/criTRia-curations.tsv). `paper/supp4_criTRia_curations.tsv` is a snapshot of that file with an added `Locus_ID` column (the STRchive locus ID), refreshed every time the pipeline runs.

## Paper items

Every file in `paper/` goes into the manuscript, and nothing else does.

| Paper item | File | Made by |
| :--- | :--- | :--- |
| Figure 1 (evidence matrix) | Made by hand: [Google Sheet](https://docs.google.com/spreadsheets/d/1VEuZqvwtQWzVSBc7Aj4FUPQABwyNTCcb8XDNCZ-xoKM/edit?usp=sharing) | – |
| Figure 2 (criTRia vs GenCC heatmap) | `paper/fig2_heatmap.pdf`, `paper/fig2_heatmap.png` | `scripts/03_figures.R` |
| Figure 3 (minimum evidence for a new TR locus) | Made by hand | – |
| Supplemental File 3 (criTRia and GenCC classifications) | `paper/supp3_dataset.tsv` | `scripts/01_download_and_match.py` |
| Supplemental File 4 (criTRia curations) | `paper/supp4_criTRia_curations.tsv` | `scripts/01_download_and_match.py` |
| Supplemental File X (discordant loci) | `paper/suppX_discordant_loci.tsv` | `scripts/02_discordance.R` |

Supplemental Files 1 (SOP) and 2 (LLM prompt) are not kept in this repository.

## Running

Rebuild everything, in order, with:

```bash
./run_all.sh
```

Requires Python 3 and R with tidyverse, ggnewscale, cowplot, gt, UpSetR and data.table.

## Layout

```
run_all.sh        runs the scripts below in order
scripts/
  01_download_and_match.py   download GenCC, STRchive and criTRia data; match GenCC diseases to loci
  02_discordance.R           discordance matrix and discordant loci table
  03_figures.R               Figure 2 and exploratory figures
  04_sensitivity.R           score-perturbation sensitivity analysis (exploratory)
  discordance.R              discordance flags, shared by 02 and 03
data/
  raw/          downloaded inputs (dated GenCC snapshot)
  processed/    intermediate outputs: gene_disease_matches.tsv, discordance_matrix.tsv
paper/          outputs used in the manuscript (see Paper items)
exploratory/    outputs not used in the manuscript
```

## Data and methods

### Sources

- [GenCC submissions](https://thegencc.org/download) (ClinGen, Ambry, G2P, PanelApp (Genomics England and Australia combined), Illumina, Labcorp, Lab MM, Myriad, Orphanet), downloaded on each run. If the download fails, the dated snapshot `data/raw/gencc-submissions_2026-10-02.tsv` is used instead.
- Locus disease IDs from the [STRchive](https://github.com/dashnowlab/STRchive) loci file. Provisional loci are excluded.
- criTRia curations from STRchive (saved to `paper/supp4_criTRia_curations.tsv`).

### Matching GenCC classifications to loci

A GenCC classification is assigned to a locus only when its disease ID (MONDO, OMIM or Orphanet) is one of the IDs listed for that locus in STRchive, or one accepted after manual review (`EXTRA_DISEASE_IDS` in `scripts/01_download_and_match.py`). If a group has several matching records, the most recent curation is kept; for records curated on the same date, Refuted or Disputed takes precedence, then the strongest classification. `Refuted` and `Disputed` are reported as `Contradictory`.

`paper/supp3_dataset.tsv` has one row per locus and group: `Locus_ID`, `Group`, `categorical_score`, and for GenCC rows the `Submitter`, matched `Disease_ID` and `Disease_Name`, and `Curation_Date` of the record used. Every GenCC record for each locus gene is listed in `data/processed/gene_disease_matches.tsv` with `status` set to `match`, `excluded` (listed in `EXCLUDED_DISEASE_IDS`) or `disease_mismatch`, so mismatches can be reviewed.

### Discordance

A locus is discordant when one classification is high (Definitive, Strong or Moderate) and another is low (Limited, Contradictory or No Known); Supportive counts as neither. "vs criTRia" means criTRia is on the opposite side from at least one GenCC group, and "within GenCC" means GenCC groups disagree among themselves. Both flags are defined in `scripts/discordance.R` and used for Figure 2 and Supplemental File X. The scenario and rationale text in `scripts/02_discordance.R` is written by hand and should be checked when classifications change.

### Manually reviewed comparisons not matched

GenCC classifications for a locus gene that were reviewed by hand and deliberately not matched to the locus disease. These already fail the disease ID match, so they need no special handling in the script.

| Locus | GenCC disease | Groups | Reason not matched |
| :--- | :--- | :--- | :--- |
| FRA2A_AFF3 | KINSSHIP syndrome (MONDO:0851095, OMIM:619297) | G2P, Labcorp, PanelApp (all Strong) | KINSSHIP is caused by de novo AFF3 degron variants that increase AFF3 and has a distinct multisystem phenotype (horseshoe kidney, mesomelic dysplasia, seizures, hypertrichosis, pulmonary involvement; PMID 33961779), whereas FRA2A silences AFF3 and causes intellectual disability. Labcorp's record also cites the FRA2A expansion paper (PMID 24763282) but is excluded on its disease label. |
| SD5_HOXD13 | Brachydactyly-syndactyly syndrome (MONDO:0012544, OMIM:610713) | G2P (Definitive), Lab MM (Strong), Orphanet (Supportive) | Caused by a contraction of the same HOXD13 polyalanine tract (15 to 8 alanines; PMID 17236141), but a different disease from the synpolydactyly type 1 caused by expansions, which criTRia curated. Treated as a separate locus-disease relationship, as for FMR1 (FXS vs FXTAS/POF1). G2P's record also cites synpolydactyly expansion papers (e.g. PMID 8817328) but is excluded on its disease label. |

### GenCC classifications based only on non-tandem-repeat variants

Some GenCC classifications that match a locus disease are based only on non-tandem-repeat (non-TR) variants in the gene. These are included in the dataset and figures like any other classification, so their scores reflect gene-level rather than repeat-level evidence. They were identified from each record's curation date (compared with the year STRchive gives for the TR being reported) and the papers it cites. Records that cite nothing, or whose citations may include repeat alleles, are not listed.

| Locus | Group | Score | Reason |
| :--- | :--- | :--- | :--- |
| CCD_RUNX2 | Labcorp | Strong | Cites only point-mutation and deletion papers |
| CPEO_POLG | Labcorp | Strong | Cites only point-mutation and Parkinson disease papers |
| CPEO_POLG | PanelApp | Strong | Cites only a review of POLG point mutations |
| DBQD2_XYLT1 | Orphanet | Supportive | Cites only a 2014 point-mutation paper; TR reported 2019 |
| DMD_DMD | ClinGen | Definitive | Cites only non-TR variant papers |
| DMD_DMD | G2P | Definitive | Curated 2015, before the TR was reported (2016) |
| DMD_DMD | Labcorp | Strong | Cites only non-TR variant papers |
| DMD_DMD | PanelApp | Strong | Cites only the GeneReviews dystrophinopathies chapter |
| FRA2A_AFF3 | G2P | Moderate | Cites only a paper on AFF3 loss-of-function and missense variants |
| EPM1_CSTB | Labcorp | Strong | Cites only point-mutation, null-allele and mouse model papers |
| HPE5_ZIC2 | G2P | Definitive | Cites only non-TR variant papers |
| HPE5_ZIC2 | Labcorp | Strong | Cites only non-TR variant papers |
| NME_NAXE | ClinGen | Definitive | Cites only point-mutation papers; TR reported 2024 |
| NME_NAXE | G2P | Strong | Curated 2017, before the TR was reported (2024) |
| NME_NAXE | PanelApp | Strong | Cites only point-mutation papers |
| SCA27B_FGF14 | Ambry | Moderate | SCA27A curation from 2018, before the TR was reported (2023) |
| SCA27B_FGF14 | Orphanet | Supportive | SCA27A curation from 2021, before the TR was reported (2023) |
| SCA27B_FGF14 | PanelApp | Strong | SCA27A curation citing only FGF14 coding-variant papers |
| SCA4_ZFHX3 | PanelApp | Strong | Cites only ZFHX3 loss-of-function neurodevelopmental disorder papers |

## Repository contributors

- **Principal Investigator:** Harriet Dashnow
- **First Author:** Macayla Ann Weiner
