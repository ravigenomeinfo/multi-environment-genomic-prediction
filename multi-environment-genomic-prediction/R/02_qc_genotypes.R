# Step 2: genotype QC (genotypes only, so no phenotype leakage) and data assembly.
if (!file.exists("config.R")) stop("Run scripts from the repository root (the folder containing config.R).", call. = FALSE)
source("config.R"); source("R/functions.R")

X_raw <- as.matrix(read.csv("data/simulated/genotypes_raw.csv", row.names = 1, check.names = FALSE))
qc <- qc_genotypes(X_raw, config$maf_min, config$snp_miss_max, config$ind_miss_max)

ph <- read.csv("data/simulated/phenotypes_long.csv")
Y <- tapply(ph$y, list(ph$line, ph$env), mean)[rownames(qc$X), ]
gt <- read.csv("data/simulated/truth_genetic_values.csv", row.names = 1)
Gtrue <- as.matrix(gt)[rownames(qc$X), colnames(Y)]
stopifnot(!anyNA(qc$X), !anyNA(Y), identical(colnames(Y), colnames(Gtrue)))

saveRDS(qc$X, "data/processed/X_qc.rds"); saveRDS(Y, "data/processed/Y.rds")
saveRDS(Gtrue, "data/processed/Gtrue.rds")
write.csv(qc$stats, "results/tables/qc_marker_stats.csv", row.names = FALSE)
write.csv(qc$summary, "results/tables/qc_summary.csv", row.names = FALSE)
print(qc$summary)
