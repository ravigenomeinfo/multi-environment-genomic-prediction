# Step 1: simulate genotypes and a multi-environment trial (MET).
# ALL DATA ARE SIMULATED. Ground truth is written to truth_*.csv files.
if (!file.exists("config.R")) stop("Run scripts from the repository root (the folder containing config.R).", call. = FALSE)
source("config.R"); source("R/functions.R")
set.seed(config$seed)

geno <- simulate_genotypes(config$n_lines, config$n_snps, config$n_chr, config$n_pops,
                           config$fst, config$missing_rate, config$ld_rho)
met <- simulate_met(geno$X_true, config)

write.csv(geno$X, "data/simulated/genotypes_raw.csv", quote = FALSE)
write.csv(transform(make_long(met$y)[, c("line", "env", "y")], y = round(y, 4)),
          "data/simulated/phenotypes_long.csv", row.names = FALSE, quote = FALSE)
write.csv(data.frame(env = rownames(met$C), round(met$C, 4)),
          "data/simulated/env_covariates.csv", row.names = FALSE, quote = FALSE)
write.csv(data.frame(line = rownames(met$g), round(met$g, 4)),
          "data/simulated/truth_genetic_values.csv", row.names = FALSE, quote = FALSE)
write.csv(data.frame(env = names(met$h2), h2_plot_level = round(met$h2, 3), env_mean = round(met$mu, 3)),
          "data/simulated/truth_env.csv", row.names = FALSE, quote = FALSE)
write.csv(data.frame(line = names(geno$pop), subpop = geno$pop),
          "data/simulated/truth_lines.csv", row.names = FALSE, quote = FALSE)

cat("Realised genetic correlations between environments (simulation truth):\n")
print(round(cor(met$g), 2))
cat("Phenotypic correlations between environments:\n"); print(round(cor(met$y), 2))
