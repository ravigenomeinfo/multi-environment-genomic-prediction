# Step 3: genomic relationship matrix G (lines) and environmental kernel K_E (environments).
if (!file.exists("config.R")) stop("Run scripts from the repository root (the folder containing config.R).", call. = FALSE)
source("config.R"); source("R/functions.R")

X <- readRDS("data/processed/X_qc.rds")
G <- vanraden_G(X)
cov_df <- read.csv("data/simulated/env_covariates.csv", row.names = 1)
KE <- env_kernel(as.matrix(cov_df))            # built from ALL 5 measured covariates
saveRDS(G, "data/processed/G.rds"); saveRDS(KE, "data/processed/KE.rds")

Rg <- cor(readRDS("data/processed/Gtrue.rds"))
off <- upper.tri(Rg)
cat("Environmental kernel (K_E):\n"); print(round(KE, 2))
cat("Correlation between K_E and realised genetic correlation (off-diagonal):",
    round(cor(KE[off], Rg[off]), 2), "\n")
cat("(Informative but imperfect: 2 of the 5 covariates are unrelated noise.)\n")
