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

echo "== 4. Sensitivity analysis (±1 point score perturbation)"
Rscript scripts/04_sensitivity.R

echo "== 5. Inter-curator agreement"
# Independent curator (30-locus re-curation sheet) vs the original criTRia curations in STRchive
mkdir -p data/curators
curl -sSfL https://raw.githubusercontent.com/dashnowlab/STRchive/refs/heads/main/data/criTRia-curations.json \
  -o data/curators/STRchive-criTRia-curations.json
python3 scripts/import_curations.py --sheet 1DF2NPHb-CHq2XrvoyxNU_5yAuHo_szQr7SubEEfddMo \
  -o data/curators/curations-gabriel.json
python3 scripts/compare_curators.py --rater gabriel data/curators/curations-gabriel.json \
  --rater group data/curators/STRchive-criTRia-curations.json -o data/curators/curator-comparison.tsv
Rscript scripts/05_curator_agreement.R

echo "Done."
