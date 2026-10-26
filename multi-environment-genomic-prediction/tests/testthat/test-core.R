small_cfg <- list(n_env = 3, n_env_cov = 4, n_env_cov_causal = 2, env_cor_bandwidth = 1.1,
                  h2_range = c(0.4, 0.6), env_mean_sd = 1, n_qtl = 30)
make_small <- function(seed = 1, n = 80) {
  set.seed(seed)
  g <- simulate_genotypes(n, 400, 4, 2, 0.1, 0.02, 0.95)
  m <- simulate_met(g$X_true, small_cfg)
  qc <- qc_genotypes(g$X, 0.05, 0.1, 0.2)
  G <- vanraden_G(qc$X)
  list(g = g, m = m, G = G, d = make_long(m$y[rownames(G), ]), KE = env_kernel(m$C))
}

test_that("simulation is reproducible and well-formed", {
  a <- make_small(3); b <- make_small(3)
  expect_identical(a$m$y, b$m$y)
  expect_equal(dim(a$m$g), c(80, 3))
  expect_equal(unname(apply(a$m$g, 2, sd)), rep(1, 3), tolerance = 1e-6)
})

test_that("environmental kernel is a valid similarity matrix", {
  s <- make_small()
  expect_equal(s$KE, t(s$KE)); expect_equal(unname(diag(s$KE)), rep(1, 3))
  expect_gte(min(eigen(s$KE, symmetric = TRUE, only.values = TRUE)$values), -1e-8)
})

test_that("long-format row order matches kronecker(K_E, G) level order", {
  s <- make_small()
  expect_identical(s$d$envline, rownames(kron_G(s$KE, s$G)))
})

test_that("M1 and M2 predictions are finite and informative for new lines (CV1)", {
  s <- make_small(5, n = 100); n <- nrow(s$G); E <- 3
  test_lines <- rownames(s$G)[1:20]
  tr <- which(!(s$d$line %in% test_lines))
  for (m in c("M1", "M2")) {
    P <- predict_met(m, s$d, tr, s$G)
    expect_equal(dim(P), c(n, E)); expect_true(all(is.finite(P)))
    r <- mean(sapply(1:E, function(e) cor(s$m$g[test_lines, e], P[test_lines, e])))
    expect_gt(r, 0)
  }
})

test_that("a completely unseen environment can be predicted by M2 and M5 (CV0)", {
  s <- make_small(7, n = 60)
  KG <- kron_G(s$KE, s$G); tr <- which(s$d$env != "E3")
  for (m in c("M2", "M5")) {
    P <- predict_met(m, s$d, tr, s$G, KG)
    expect_true(all(is.finite(P[, "E3"])))
  }
})

test_that("cell_scores skips environments with too few test cells", {
  Y <- matrix(rnorm(60), 20, 3, dimnames = list(NULL, c("E1", "E2", "E3")))
  mask <- matrix(FALSE, 20, 3, dimnames = dimnames(Y)); mask[1:10, "E1"] <- TRUE; mask[1:2, "E2"] <- TRUE
  out <- cell_scores(Y, Y, Y, mask, data.frame(design = "x"))
  expect_equal(out$env, "E1"); expect_equal(out$r, 1)
})
