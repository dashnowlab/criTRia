# criTRia
criTRia (pronounce criteria) is a formalized scoring framework designed to accurately evaluate TR-disease relationships. criTRia builds on current best-practice scoring procedures developed by the Clinical Genome Resource (ClinGen) while (re)curating based on locus-level (rather than gene-level) classification, introducing TR-specific evidence categories and reweighted scoring. In this repository you will find all data and code used to evaluate success and create the figures for the Mnauscript and SOP, up to date as of 4/22/26.

## Download and format Google Sheet (Python)

Use `download_sheet_to_criteria_dataset.py` to download a Google Sheet tab and standardize it to the same format as `criTRia_Dataset.csv`:

- Output columns: `Gene,Group,categorical_score`
- `Refuted` and `Disputed` scores are converted to `Contradictory`

Run from the repository root:

```bash
python3 download_sheet_to_criteria_dataset.py
```

Write directly to a specific file (for example, replacing the dataset):

```bash
python3 download_sheet_to_criteria_dataset.py --output criTRia_Dataset.csv
```

Use a different Google Sheet URL:

```bash
python3 download_sheet_to_criteria_dataset.py --sheet-url "<google-sheet-url>"
```
