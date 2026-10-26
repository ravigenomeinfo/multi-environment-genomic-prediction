# Shared functions. Genotypes coded -1/0/1; simulated inbred lines are -1 or 1.
suppressPackageStartupMessages({ library(rrBLUP); library(sommer); library(ggplot2) })

# ------------------------------------------------------------- simulation
simulate_genotypes <- function(n_lines, n_snps, n_chr, n_pops, fst, missing_rate, ld_rho) {
  # Balding-Nichols subpopulations + Gaussian-copula LD along chromosomes
  # (simplified; not a coalescent simulation).
  p_anc <- pmin(pmax(rbeta(n_snps, 0.7, 0.7), 0.02), 0.98)
  a <- p_anc * (1 - fst) / fst; b <- (1 - p_anc) * (1 - fst) / fst
  p_pop <- sapply(seq_len(n_pops), function(k) rbeta(n_snps, a, b))
  pop <- sample(rep_len(seq_len(n_pops), n_lines))
  chr <- sort(rep_len(seq_len(n_chr), n_snps))
  Z <- matrix(rnorm(n_lines * n_snps), n_lines, n_snps)
  for (j in seq_len(n_snps)[-1]) if (chr[j] == chr[j - 1]) Z[, j] <- ld_rho * Z[, j - 1] + sqrt(1 - ld_rho^2) * Z[, j]
  X <- ifelse(Z < t(qnorm(p_pop[, pop, drop = FALSE])), 1, -1)
  dimnames(X) <- list(sprintf("L%03d", seq_len(n_lines)), sprintf("SNP%04d", seq_len(n_snps)))
  X_true <- X
  X[sample(length(X), round(missing_rate * length(X)))] <- NA
  list(X = X, X_true = X_true, pop = setNames(pop, rownames(X)))
}

# Environmental covariates -> Gaussian kernel (similarity between environments)
env_kernel <- function(C, bandwidth = NULL) {
  D2 <- as.matrix(dist(scale(C)))^2
  if (is.null(bandwidth)) bandwidth <- median(D2[upper.tri(D2)])
  K <- exp(-D2 / bandwidth)
  dimnames(K) <- list(rownames(C), rownames(C)); K
}

# Genetic values per environment with known between-environment genetic correlation.
# Correlation is generated from the first `n_causal` covariates only, so a kernel built
# from ALL measured covariates is informative but imperfect (noise covariates included).
simulate_met <- function(X_true, cfg) {
  E <- cfg$n_env; env <- paste0("E", seq_len(E))
  C <- matrix(rnorm(E * cfg$n_env_cov), E, cfg$n_env_cov, dimnames = list(env, paste0("cov", seq_len(cfg$n_env_cov))))
  Ctrue <- C[, seq_len(cfg$n_env_cov_causal), drop = FALSE]
  Rg <- exp(-(as.matrix(dist(Ctrue))^2) / (2 * cfg$env_cor_bandwidth^2 * cfg$n_env_cov_causal))
  dimnames(Rg) <- list(env, env)

  p <- (colMeans(X_true) + 1) / 2; maf <- pmin(p, 1 - p)
  qtl <- sample(which(maf > 0.10), cfg$n_qtl)
  B <- matrix(rnorm(cfg$n_qtl * E), cfg$n_qtl, E) %*% chol(Rg)  # effects correlated across env
  g <- scale(X_true[, qtl, drop = FALSE] %*% B)                  # unit genetic variance per env
  h2 <- seq(cfg$h2_range[1], cfg$h2_range[2], length.out = E)
  mu <- rnorm(E, 8, cfg$env_mean_sd)
  sd_e <- sqrt((1 - h2) / h2)
  y <- sapply(seq_len(E), function(e) mu[e] + g[, e] + rnorm(nrow(g), 0, sd_e[e]))
  g <- matrix(g, nrow(g), E)
  dimnames(g) <- dimnames(y) <- list(rownames(X_true), env)
  list(y = y, g = g, C = C, Rg_design = Rg, h2 = setNames(h2, env), mu = setNames(mu, env), qtl = qtl)
}

# --------------------------------------------------------------------- QC
qc_genotypes <- function(X, maf_min, snp_miss_max, ind_miss_max) {
  n0 <- nrow(X); m0 <- ncol(X)
  X <- X[rowMeans(is.na(X)) <= ind_miss_max, , drop = FALSE]
  p <- (colMeans(X, na.rm = TRUE) + 1) / 2
  stats <- data.frame(snp = colnames(X), maf = pmin(p, 1 - p), missing = colMeans(is.na(X)))
  stats$kept <- stats$missing <= snp_miss_max & stats$maf >= maf_min
  X <- X[, stats$kept, drop = FALSE]
  n_imp <- sum(is.na(X)); cm <- colMeans(X, na.rm = TRUE)
  for (j in seq_len(ncol(X))) X[is.na(X[, j]), j] <- cm[j]
  list(X = X, stats = stats,
       summary = data.frame(step = c("lines_before", "lines_after", "snps_before", "snps_after", "calls_imputed"),
                            value = c(n0, nrow(X), m0, ncol(X), n_imp)))
}

vanraden_G <- function(X) {
  p <- (colMeans(X) + 1) / 2
  Z <- sweep(X, 2, 2 * p - 1)
  G <- tcrossprod(Z) / (2 * sum(p * (1 - p)))
  dimnames(G) <- list(rownames(X), rownames(X)); G
}

# ------------------------------------------------------- long-format data
# Row order is environment-major (all lines in E1, then E2, ...), matching kronecker(K_E, G).
make_long <- function(y_mat) {
  d <- data.frame(line = rep(rownames(y_mat), times = ncol(y_mat)),
                  env  = rep(colnames(y_mat), each = nrow(y_mat)),
                  y    = as.vector(y_mat), stringsAsFactors = FALSE)
  d$envline <- paste(d$env, d$line, sep = ":")
  d
}

kron_G <- function(KE, G) {
  KG <- kronecker(KE, G)
  lv <- paste(rep(rownames(KE), each = nrow(G)), rep(rownames(G), times = nrow(KE)), sep = ":")
  dimnames(KG) <- list(lv, lv); KG
}

# ----------------------------------------------------------------- models
# M1 single-environment GBLUP (one model per environment, no borrowing of information)
# M2 across-environment, genotype main effect only
# M3 main + environment-specific genomic variances (diagonal G x E)
# M4 unstructured genetic covariance across environments
# M5 main + G x E with an environmental-covariate kernel (K_E (x) G)
model_labels <- c(M1 = "M1 single-env", M2 = "M2 main effect", M3 = "M3 main + diag GxE",
                  M4 = "M4 unstructured", M5 = "M5 main + kernel GxE")

# Heterogeneous residual variance per environment. If the inversion tolerance reports a
# singular system, refit with progressively looser tolerances (1e-3, 1e-2, 1e-1); the
# unstructured model starts at 1e-3 because it needs it on small panels.
# sommer evaluates Gu = G / KG by NAME in the calling scope, so the matrices are passed
# in as arguments named G and KG (not found via globals).
fit_mmes <- function(rhs_random, d_fit, tol = 1e-6, G = NULL, KG = NULL) {
  run <- function(tol) mmes(y ~ env, random = rhs_random, rcov = ~ vsm(dsm(env), ism(units)),
                            data = d_fit, tolParInv = tol, dateWarning = FALSE, verbose = FALSE)
  for (t in sort(unique(c(tol, 1e-3, 1e-2, 1e-1)[c(tol, 1e-3, 1e-2, 1e-1) >= tol]))) {
    fit <- tryCatch(run(t), error = function(e) e)
    if (!inherits(fit, "error")) return(fit)
    if (!grepl("singular", conditionMessage(fit))) stop(fit)
  }
  stop(fit)
}

# Fixed part per environment. Training environments use their estimate; an unseen
# environment gets the mean over training environments (this affects bias, not the
# within-environment correlations reported here).
env_fixed <- function(fit, envs_all, envs_train) {
  b <- drop(fit$b); names(b) <- rownames(fit$b)
  eff <- setNames(rep(0, length(envs_train)), envs_train)
  for (e in setdiff(envs_train, envs_train[1])) eff[e] <- b[e]
  base <- b["Intercept"] + eff
  out <- setNames(rep(mean(base), length(envs_all)), envs_all)
  out[envs_train] <- base
  out
}

# Fit on rows `train_idx` of d; predict every (line, environment) cell.
# Returns a lines x environments matrix.
predict_met <- function(model, d, train_idx, G, KG = NULL) {
  lines <- rownames(G); envs <- sort(unique(d$env)); n <- length(lines); E <- length(envs)
  tr <- d[train_idx, ]
  tr$env <- factor(tr$env, levels = envs[envs %in% tr$env])
  tr$line <- factor(tr$line, levels = lines)
  envs_tr <- levels(tr$env)
  P <- matrix(NA_real_, n, E, dimnames = list(lines, envs))

  if (model == "M1") {
    for (e in envs_tr) {
      te <- tr[tr$env == e, ]; obs <- as.character(te$line)
      fit <- mixed.solve(te$y, K = G[obs, obs])
      # lines unobserved in this environment: BLUP through their covariance with observed lines
      u_all <- G[lines, obs, drop = FALSE] %*% solve(G[obs, obs] + diag(1e-6, length(obs)), fit$u)
      P[, e] <- as.numeric(fit$beta) + as.numeric(u_all)
    }
    return(P)
  }
  if (model == "M2") {
    fit <- fit_mmes(~ vsm(ism(line), Gu = G), tr, G = G, KG = KG)
    u <- fit$uList[[1]][lines, 1]; fx <- env_fixed(fit, envs, envs_tr)
    for (e in envs) P[, e] <- fx[e] + u
  } else if (model == "M3") {
    fit <- fit_mmes(~ vsm(ism(line), Gu = G) + vsm(dsm(env), ism(line), Gu = G), tr, G = G, KG = KG)
    u <- fit$uList[[1]][lines, 1]; ug <- fit$uList[[2]]; fx <- env_fixed(fit, envs, envs_tr)
    for (e in envs) P[, e] <- fx[e] + u + (if (e %in% colnames(ug)) ug[lines, e] else 0)
  } else if (model == "M4") {
    fit <- fit_mmes(~ vsm(usm(env), ism(line), Gu = G), tr, tol = 1e-3, G = G)
    ug <- fit$uList[[1]]; fx <- env_fixed(fit, envs, envs_tr)
    for (e in envs) P[, e] <- fx[e] + (if (e %in% colnames(ug)) ug[lines, e] else rowMeans(ug[lines, envs_tr, drop = FALSE]))
  } else if (model == "M5") {
    tr$envline <- factor(tr$envline, levels = rownames(KG))
    fit <- fit_mmes(~ vsm(ism(line), Gu = G) + vsm(ism(envline), Gu = KG), tr, G = G, KG = KG)
    u <- fit$uList[[1]][lines, 1]; uk <- fit$uList[[2]][, 1]; fx <- env_fixed(fit, envs, envs_tr)
    for (k in seq_along(envs)) P[, envs[k]] <- fx[envs[k]] + u + uk[((k - 1) * n + 1):(k * n)]
  } else stop("unknown model: ", model)
  P
}

# ----------------------------------------------------------------- metrics
# Within-environment correlation over the test cells of each environment.
cell_scores <- function(P, Ymat, Gtrue, test_mask, meta) {
  do.call(rbind, lapply(colnames(P), function(e) {
    i <- which(test_mask[, e])
    if (length(i) < 5) return(NULL)
    data.frame(meta, env = e, n_test = length(i),
               r = suppressWarnings(cor(Ymat[i, e], P[i, e])),
               r_tgv = suppressWarnings(cor(Gtrue[i, e], P[i, e])))
  }))
}

# ------------------------------------------------------------------- plots
okabe_ito <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
theme_paper <- function() theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"))
save_fig <- function(plot, name, width = 6, height = 4.5)
  ggsave(file.path("results", "figures", name), plot, width = width, height = height, dpi = 300)
check_repo_root <- function() if (!file.exists("config.R") || !dir.exists("R"))
  stop("Run scripts from the repository root (the folder containing config.R).", call. = FALSE)
