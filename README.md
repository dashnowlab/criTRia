# criTRia

criTRia (pronounce criteria) is a formalized scoring framework designed to accurately evaluate TR-disease relationships. criTRia builds on current best-practice scoring procedures developed by the Clinical Genome Resource (ClinGen) while (re)curating based on locus-level (rather than gene-level) classification, introducing TR-specific evidence categories and reweighted scoring. In this repository you will find all data and code used to evaluate success and create the figures for the manuscript and SOP, up to date as of the latest commit.

## Download and format Google Sheet (Python)

Use `download_sheet_to_criteria_dataset.py` to download a Google Sheet tab and standardize it to the same format as `criTRia_Dataset.csv`:

- Output columns: `Gene,Group,categorical_score`
- `Refuted` and `Disputed` scores are converted to `Contradictory`

Run from the repository root:

```bash
python3 download_sheet_to_criteria_dataset.py
```

Merge data and generate figures

```bash
Rscript criTRia_Figure_Script.R
```

## Figure source files

This repository includes links to editable source files for figures used in the criTRia manuscript.

| Figure # | Description | Editable Source Link |
| :--- | :--- | :--- |
| **Figure 1** | criTRia scoring framework overview | [https://docs.google.com/spreadsheets/d/1VEuZqvwtQWzVSBc7Aj4FUPQABwyNTCcb8XDNCZ-xoKM/edit?usp=sharing](https://docs.google.com/spreadsheets/d/1VEuZqvwtQWzVSBc7Aj4FUPQABwyNTCcb8XDNCZ-xoKM/edit?usp=sharing) |
| **Figure 2** | Comparison of criTRia vs Gene Curation Coalition | [https://docs.google.com/spreadsheets/d/1VEuZqvwtQWzVSBc7Aj4FUPQABwyNTCcb8XDNCZ-xoKM/edit?usp=sharing](https://docs.google.com/spreadsheets/d/1VEuZqvwtQWzVSBc7Aj4FUPQABwyNTCcb8XDNCZ-xoKM/edit?usp=sharing) |
| **Figure 3** | TR-disease association results | [https://github.com/dashnowlab/criTRia/blob/main/criTRia%20Figure%20Script.R](https://github.com/dashnowlab/criTRia/blob/main/criTRia%20Figure%20Script.R) |

## Contact information

- **Principal Investigator:** Harriet Dashnow
- **Author:** Macayla Ann Weiner
