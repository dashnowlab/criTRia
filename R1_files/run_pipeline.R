# run_pipeline.R
# Master orchestration script for the criTRia disease-matched re-curation pipeline.
# Executes all steps completely in R.

cat("=======================================================================\n")
cat("Starting criTRia Disease-Matched Re-Curation & Analysis Pipeline (R Only)\n")
cat("=======================================================================\n\n")

# Step 1: Audit and Match GenCC submissions
cat(">> Executing Step 1: 01_audit_and_match_gencc.R\n")
source("01_audit_and_match_gencc.R")
cat("\n-----------------------------------------------------------------------\n\n")

# Step 2: Build Discordance Matrix and Supplemental Table S3b
cat(">> Executing Step 2: 02_build_discordance_matrix.R\n")
source("02_build_discordance_matrix.R")
cat("\n-----------------------------------------------------------------------\n\n")

# Step 3: Plot publication-ready Figure 2 Heatmap
cat(">> Executing Step 3: 03_plot_figure2_matched.R\n")
source("03_plot_figure2_matched.R")
cat("\n-----------------------------------------------------------------------\n\n")

# Step 4: Generate Supplemental File 6 R2 (Discordant Loci Breakdown)
cat(">> Executing Step 4: 06_build_discordant_table.R\n")
source("06_build_discordant_table.R")
cat("\n-----------------------------------------------------------------------\n\n")

# Step 5: Generate Supplemental File 3 R2 (Clean Long-Format Dataset)
cat(">> Executing Step 5: 07_generate_supplemental_file_3.R\n")
source("07_generate_supplemental_file_3.R")
cat("\n-----------------------------------------------------------------------\n\n")

# Step 6: Run Score Perturbation Sensitivity Analysis
cat(">> Executing Step 6: 08_run_sensitivity_analysis.R\n")
source("08_run_sensitivity_analysis.R")
cat("\n=======================================================================\n")
cat("Pipeline execution complete! All deliverables generated successfully.\n")
cat("=======================================================================\n")
