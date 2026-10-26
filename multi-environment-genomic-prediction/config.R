# Central configuration. Every script sources this file.
config <- list(
  seed = 2026,

  # ---- simulated multi-environment trial (MET) ----
  n_lines      = 120,
  n_snps       = 1000,
  n_chr        = 10,
  n_pops       = 3,
  fst          = 0.10,
  ld_rho       = 0.97,
  missing_rate = 0.02,
  n_qtl        = 60,
  n_env        = 5,
  n_env_cov    = 5,      # measured environmental covariates per environment
  n_env_cov_causal = 3,  # only the first 3 drive genetic correlations between environments
  env_cor_bandwidth = 1.1,  # larger = more similar environments (higher genetic correlations)
  h2_range     = c(0.35, 0.65),  # per-environment heritability (plot-level)
  env_mean_sd  = 1.0,

  # ---- genotype QC ----
  maf_min = 0.05, snp_miss_max = 0.10, ind_miss_max = 0.20,

  # ---- validation ----
  k_folds = 5,
  n_reps  = 2
)

# QUICK=1 Rscript run_all.R  -> one repetition (about half the run time)
if (nzchar(Sys.getenv("QUICK"))) config$n_reps <- 1
