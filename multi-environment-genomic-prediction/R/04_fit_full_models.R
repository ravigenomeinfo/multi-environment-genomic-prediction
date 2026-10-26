# Step 4: fit the multi-environment models on ALL data (descriptive, not validation).
# Compares model fit (AIC/BIC) and checks whether the unstructured model (M4)
# recovers the simulated genetic correlations between environments.
if (!file.exists("config.R")) stop("Run scripts from the repository root (the folder containing config.R).", call. = FALSE)
source("config.R"); source("R/functions.R")

G <- readRDS("data/processed/G.rds"); KE <- readRDS("data/processed/KE.rds")
Y <- readRDS("data/processed/Y.rds"); Gtrue <- readRDS("data/processed/Gtrue.rds")
d <- make_long(Y); d$env <- factor(d$env); d$line <- factor(d$line, levels = rownames(G))
KG <- kron_G(KE, G); d$envline <- factor(d$envline, levels = rownames(KG))

fits <- list(
  M2 = fit_mmes(~ vsm(ism(line), Gu = G), d, G = G, KG = KG),
  M3 = fit_mmes(~ vsm(ism(line), Gu = G) + vsm(dsm(env), ism(line), Gu = G), d, G = G, KG = KG),
  M4 = fit_mmes(~ vsm(usm(env), ism(line), Gu = G), d, tol = 1e-3, G = G),
  M5 = fit_mmes(~ vsm(ism(line), Gu = G) + vsm(ism(envline), Gu = KG), d, G = G, KG = KG))

fit_tab <- do.call(rbind, lapply(names(fits), function(m)
  data.frame(model = m, description = model_labels[[m]],
             logLik = round(tail(as.numeric(fits[[m]]$llik), 1), 2),
             AIC = round(as.numeric(fits[[m]]$AIC)[1], 2), BIC = round(as.numeric(fits[[m]]$BIC)[1], 2),
             converged = isTRUE(fits[[m]]$convergence))))
write.csv(fit_tab, "results/tables/model_fit_full_data.csv", row.names = FALSE)
cat("Model fit on all data (lower AIC/BIC = better; includes in-sample fit only):\n")
print(fit_tab)

# Unstructured genetic (co)variance -> correlation matrix
vc <- summary(fits$M4)$varcomp
envs <- colnames(Y); S <- matrix(0, length(envs), length(envs), dimnames = list(envs, envs))
for (nm in rownames(vc)) {
  p <- strsplit(nm, ":")[[1]]
  if (length(p) >= 5 && p[1] == "env" && p[2] == "line") { S[p[4], p[5]] <- vc[nm, "VarComp"]; S[p[5], p[4]] <- vc[nm, "VarComp"] }
}
R_est <- suppressWarnings(cov2cor(S)); R_true <- cor(Gtrue)
write.csv(round(S, 4), "results/tables/M4_genetic_covariance.csv")
write.csv(round(R_est, 3), "results/tables/M4_genetic_correlation_estimated.csv")
write.csv(round(R_true, 3), "results/tables/genetic_correlation_realised_truth.csv")
off <- upper.tri(R_true)
cat("\nEstimated vs realised genetic correlation (off-diagonal): r =",
    round(cor(R_est[off], R_true[off], use = "complete.obs"), 2),
    "| mean abs error =", round(mean(abs(R_est[off] - R_true[off]), na.rm = TRUE), 2), "\n")
ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
cat("\nSmallest eigenvalue of the M4 genetic covariance matrix:", signif(min(ev), 3),
    if (min(ev) < 0) "(negative: estimate is not a valid covariance matrix; unstructured estimation is unreliable at this sample size)" else "", "\n")
write.csv(data.frame(eigenvalue = ev), "results/tables/M4_covariance_eigenvalues.csv", row.names = FALSE)
cat("\nEstimated (M4) genetic correlations:\n"); print(round(R_est, 2))
cat("Realised (truth):\n"); print(round(R_true, 2))
saveRDS(list(R_est = R_est, R_true = R_true, KE = KE), "data/processed/cor_matrices.rds")
