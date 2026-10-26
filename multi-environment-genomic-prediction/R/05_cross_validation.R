# Step 5: three validation designs (Burgueno et al. 2012 style).
#   CV1: predict NEW lines (never observed in any environment)
#   CV2: predict lines observed in other environments but not in the target one
#        (random cells masked, as in incomplete field trials)
#   CV0: predict a NEW ENVIRONMENT (leave one environment out; needs environmental information)
# Accuracy = within-environment correlation between predictions and (i) observed
# phenotypes (r) and (ii) simulated true genetic values (r_tgv; simulation only).
# Results are written after every fold, so an interrupted run keeps its progress.
if (!file.exists("config.R")) stop("Run scripts from the repository root (the folder containing config.R).", call. = FALSE)
source("config.R"); source("R/functions.R")

G <- readRDS("data/processed/G.rds"); KE <- readRDS("data/processed/KE.rds")
Y <- readRDS("data/processed/Y.rds"); Gtrue <- readRDS("data/processed/Gtrue.rds")
n <- nrow(Y); E <- ncol(Y); envs <- colnames(Y)
d <- make_long(Y); KG <- kron_G(KE, G)
stopifnot(identical(d$envline, rownames(KG)))   # row order must match kronecker(K_E, G)

out_file <- "results/tables/cv_metrics_all.csv"
unlink("results/tables/cv_failures.csv")   # remove stale file from earlier runs
results <- list(); failures <- list()
add_run <- function(design, rep, fold, model, train_idx, test_mask) {
  t0 <- Sys.time()
  P <- tryCatch(predict_met(model, d, train_idx, G, KG), error = function(e) {
    failures[[length(failures) + 1]] <<- data.frame(design, rep, fold, model, error = conditionMessage(e)); NULL })
  if (!is.null(P)) {
    results[[length(results) + 1]] <<- cell_scores(P, Y, Gtrue, test_mask, data.frame(design, rep, fold, model))
    write.csv(do.call(rbind, results), out_file, row.names = FALSE)
  }
  message(sprintf("%s rep %s fold %s %s: %.0fs%s", design, rep, fold, model,
                  as.numeric(Sys.time() - t0, units = "secs"), if (is.null(P)) "  FAILED" else ""))
}
models_all <- c("M1", "M2", "M3", "M4", "M5")

# ---- CV1: new lines ----
set.seed(config$seed + 1)
for (rep in seq_len(config$n_reps)) {
  fl <- sample(rep_len(seq_len(config$k_folds), n))
  for (f in seq_len(config$k_folds)) {
    test_lines <- rownames(Y)[fl == f]
    mask <- matrix(rownames(Y) %in% test_lines, n, E, dimnames = dimnames(Y))
    for (m in models_all) add_run("CV1_new_lines", rep, f, m, which(!(d$line %in% test_lines)), mask)
  }
}

# ---- CV2: lines observed elsewhere, missing in the target environment ----
set.seed(config$seed + 2)
for (rep in seq_len(config$n_reps)) {
  cell_fold <- matrix(sample(rep_len(seq_len(config$k_folds), n * E)), n, E, dimnames = dimnames(Y))
  for (f in seq_len(config$k_folds))
    for (m in models_all) add_run("CV2_tested_elsewhere", rep, f, m, which(as.vector(cell_fold) != f), cell_fold == f)
}

# ---- CV0: new environment ----
# M1 cannot predict an untested environment. M3/M4 have no covariance linking an unseen
# environment to the others, so they reduce to the main effect (= M2) and are omitted.
for (e in envs) {
  mask <- matrix(FALSE, n, E, dimnames = dimnames(Y)); mask[, e] <- TRUE
  for (m in c("M2", "M5")) add_run("CV0_new_environment", 1, e, m, which(d$env != e), mask)
}

res <- do.call(rbind, results)
write.csv(res, out_file, row.names = FALSE)
if (length(failures)) write.csv(do.call(rbind, failures), "results/tables/cv_failures.csv", row.names = FALSE)

summ <- aggregate(cbind(r, r_tgv) ~ design + model, data = res, FUN = function(z) mean(z, na.rm = TRUE))
summ$sd_r <- aggregate(r ~ design + model, data = res, FUN = function(z) sd(z, na.rm = TRUE))$r
fits_done <- unique(res[, c("design", "model", "rep", "fold")])
summ$n_fits <- as.vector(table(factor(paste(fits_done$design, fits_done$model), levels = paste(summ$design, summ$model))))
summ[, c("r", "r_tgv", "sd_r")] <- round(summ[, c("r", "r_tgv", "sd_r")], 3)
write.csv(summ, "results/tables/cv_summary.csv", row.names = FALSE)
by_env <- aggregate(cbind(r, r_tgv) ~ design + model + env, data = res, FUN = function(z) mean(z, na.rm = TRUE))
write.csv(by_env, "results/tables/cv_summary_by_env.csv", row.names = FALSE)
print(summ)
if (length(failures)) message("Some fits failed; see results/tables/cv_failures.csv")
