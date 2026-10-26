# Step 5b: paired model comparison. For each validation design, the accuracy of each model
# minus the accuracy of the main-effect model M2, computed on the same fold/environment
# cells, averaged over environments within each repetition x fold, then summarised as
# mean +/- standard error over repetition x fold units (CV0: over held-out environments).
if (!file.exists("config.R")) stop("Run scripts from the repository root (the folder containing config.R).", call. = FALSE)
res <- read.csv("results/tables/cv_metrics_all.csv")
wide <- reshape(res[, c("design", "rep", "fold", "env", "model", "r")],
                idvar = c("design", "rep", "fold", "env"), timevar = "model", direction = "wide")
out <- list()
for (dsg in unique(wide$design)) for (m in setdiff(sub("^r\\.", "", grep("^r\\.", names(wide), value = TRUE)), "M2")) {
  w <- wide[wide$design == dsg, ]
  col <- paste0("r.", m); if (!col %in% names(w) || all(is.na(w[[col]]))) next
  w$diff <- w[[col]] - w[["r.M2"]]
  units <- aggregate(diff ~ rep + fold, data = w, FUN = mean)       # average over environments
  if (dsg == "CV0_new_environment") units <- aggregate(diff ~ fold, data = w, FUN = mean)
  out[[length(out) + 1]] <- data.frame(design = dsg, model = m, mean_diff_vs_M2 = mean(units$diff),
    se = sd(units$diff) / sqrt(nrow(units)), n_units = nrow(units),
    frac_units_better = mean(units$diff > 0))
}
tab <- do.call(rbind, out); tab[, 3:4] <- round(tab[, 3:4], 3); tab$frac_units_better <- round(tab$frac_units_better, 2)
write.csv(tab, "results/tables/cv_paired_vs_M2.csv", row.names = FALSE)
print(tab)
