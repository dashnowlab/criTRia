# criTRia

criTRia (pronounce criteria) is a formalized scoring framework designed to accurately evaluate TR-disease relationships. criTRia builds on current best-practice scoring procedures developed by the Clinical Genome Resource (ClinGen) while (re)curating based on locus-level (rather than gene-level) classification, introducing TR-specific evidence categories and reweighted scoring. In this repository you will find all data and code used to evaluate success and create the figures for the manuscript and SOP, up to date as of the latest commit.

The current version of the criTRia curations is on the [STRchive criTRia page](https://strchive.org/critria), which lists each locus-disease classification with its score and links to the criTRia SOP. The curations are maintained in the STRchive GitHub repository as [`data/criTRia-curations.tsv`](https://github.com/dashnowlab/STRchive/blob/main/data/criTRia-curations.tsv). `criTRia-curations.tsv` in this repository is a snapshot of that file, refreshed every time the download script below is run.

## Download and format classifications (Python)

Use `download_sheet_to_criteria_dataset.py` to build `criTRia_Dataset.csv` from:

- [GenCC submissions](https://thegencc.org/download) (ClinGen, Ambry, G2P, PanelApp (Genomics England and Australia combined), Illumina, Labcorp, Lab MM, Myriad, Orphanet)
- criTRia curations and locus disease IDs from [STRchive](https://github.com/dashnowlab/STRchive) (a copy of the curations is saved to `criTRia-curations.tsv` on every run)

A GenCC classification is assigned to a locus only when its disease ID (MONDO, OMIM or Orphanet) is one of the IDs listed for that locus in STRchive. If a group has several matching records, the strongest is kept, with Refuted or Disputed taking precedence over any positive classification. Every record for each locus gene is written to `gene_disease_matches.tsv` with `status` set to `match`, `excluded` (listed in `EXCLUDED_DISEASE_IDS`) or `disease_mismatch`, so mismatches can be reviewed and, if they are the same disease, added to `EXTRA_DISEASE_IDS` in the script.

- Output columns: `Gene,Group,categorical_score`
- `Refuted` and `Disputed` scores are converted to `Contradictory`

Run from the repository root:

```bash
python3 download_sheet_to_criteria_dataset.py
```

### GenCC classifications based only on non-tandem-repeat variants

Some GenCC classifications that match a locus disease are based only on non-tandem-repeat (non-TR) variants in the gene. These are included in `criTRia_Dataset.csv` and the figures like any other classification, so their scores reflect gene-level rather than repeat-level evidence. They were identified from each record's curation date (compared with the year STRchive gives for the TR being reported) and the papers it cites. Records that cite nothing, or whose citations may include repeat alleles, are not listed.

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
| EPM1_CSTB | Labcorp | Strong | Cites only point-mutation, null-allele and mouse model papers |
| HPE5_ZIC2 | G2P | Definitive | Cites only non-TR variant papers |
| HPE5_ZIC2 | Labcorp | Strong | Cites only non-TR variant papers |
| NME_NAXE | G2P | Strong | Curated 2017, before the TR was reported (2024) |
| NME_NAXE | PanelApp | Strong | Cites only point-mutation papers |
| SCA27B_FGF14 | Ambry | Moderate | SCA27A curation from 2018, before the TR was reported (2023) |
| SCA27B_FGF14 | Orphanet | Supportive | SCA27A curation from 2021, before the TR was reported (2023) |
| SCA27B_FGF14 | PanelApp | Strong | SCA27A curation citing only FGF14 coding-variant papers |
| SCA4_ZFHX3 | PanelApp | Strong | Cites only ZFHX3 loss-of-function neurodevelopmental disorder papers |

Merge data and generate figures

```bash
Rscript criTRia_Figure_Script.R
```

## Figure source files

This repository includes links to editable source files for figures used in the criTRia manuscript.

| Figure # | Description | Editable Source Link |
| :--- | :--- | :--- |
| **Figure 1** | criTRia scoring framework overview | [https://docs.google.com/spreadsheets/d/1VEuZqvwtQWzVSBc7Aj4FUPQABwyNTCcb8XDNCZ-xoKM/edit?usp=sharing](https://docs.google.com/spreadsheets/d/1VEuZqvwtQWzVSBc7Aj4FUPQABwyNTCcb8XDNCZ-xoKM/edit?usp=sharing) |
| **Figure 2** | Comparison of criTRia vs Gene Curation Coalition | [https://docs.google.com/spreadsheets/d/1DonSiPVjeQLsB8HoFzn50jFCnWvMPxVZ124SzS0s7r0/edit?gid=0#gid=0](https://docs.google.com/spreadsheets/d/1DonSiPVjeQLsB8HoFzn50jFCnWvMPxVZ124SzS0s7r0/edit?gid=0#gid=0) |
| **Figure 3** | TR-disease association results | [https://github.com/dashnowlab/criTRia/blob/main/criTRia_Figure_Script.R](https://github.com/dashnowlab/criTRia/blob/main/criTRia_Figure_Script.R) |

## Repository contributors:

- **Principal Investigator:** Harriet Dashnow
- **First Author:** Macayla Ann Weiner
