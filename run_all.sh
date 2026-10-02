#!/usr/bin/env bash
# Rebuild every dataset, table and figure, in order. Run from anywhere:
#   ./run_all.sh
# Paper outputs are written to paper/, intermediate data to data/processed/,
# and outputs not used in the paper to exploratory/. See README.md.
set -euo pipefail
cd "$(dirname "$0")"

echo "== 1. Download GenCC, STRchive and criTRia data and match diseases"
python3 scripts/01_download_and_match.py

echo "== 2. Discordance matrix and discordant loci table"
Rscript scripts/02_discordance.R

echo "== 3. Figures"
Rscript scripts/03_figures.R

echo "== 4. Sensitivity analysis (exploratory)"
Rscript scripts/04_sensitivity.R

echo "Done."
