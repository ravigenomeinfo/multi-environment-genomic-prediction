# Full workflow from the repository root:  Rscript run_all.R
# Use  QUICK=1 Rscript run_all.R  for one CV repetition (faster).
source("R/functions.R"); check_repo_root()
for (s in c("R/01_generate_data.R", "R/02_qc_genotypes.R", "R/03_relationship_matrices.R",
            "R/04_fit_full_models.R", "R/05_cross_validation.R", "R/05b_paired_comparison.R", "R/06_visualize_results.R")) {
  message("\n=== ", s, " ===")
  source(s, echo = FALSE)
}
writeLines(capture.output(sessionInfo()), "results/session_info.txt")
message("\nDone. See results/tables and results/figures.")
