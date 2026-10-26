# Step 6: figures (300 dpi PNG).
if (!file.exists("config.R")) stop("Run scripts from the repository root (the folder containing config.R).", call. = FALSE)
source("config.R"); source("R/functions.R")

res    <- read.csv("results/tables/cv_metrics_all.csv")
by_env <- read.csv("results/tables/cv_summary_by_env.csv")
cors   <- readRDS("data/processed/cor_matrices.rds")
design_lab <- c(CV1_new_lines = "CV1: new lines",
                CV2_tested_elsewhere = "CV2: tested in other environments",
                CV0_new_environment = "CV0: new environment")
prep <- function(df) { df$design <- factor(design_lab[df$design], levels = design_lab)
                       df$model <- factor(model_labels[df$model], levels = model_labels); df }
by_env <- prep(by_env); res <- prep(res)
pal <- setNames(okabe_ito[1:5], model_labels)

# Fig 1: genetic correlations between environments (truth, M4 estimate, covariate kernel)
envs <- rownames(cors$R_true)
mk <- function(M, lab) { d <- expand.grid(e1 = envs, e2 = envs); d$value <- as.vector(M[envs, envs]); d$panel <- lab; d }
f1d <- rbind(mk(cors$R_true, "Simulated truth\n(realised genetic correlation)"),
             mk(pmax(pmin(cors$R_est, 1), -1), "M4 unstructured estimate\n(clipped to [-1, 1])"),
             mk(cors$KE, "Environmental kernel K_E\n(from covariates)"))
f1d$panel <- factor(f1d$panel, levels = unique(f1d$panel))
f1 <- ggplot(f1d, aes(e1, e2, fill = value)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
  scale_fill_gradient(low = "white", high = "#0072B2", limits = c(0, 1), na.value = "grey85", name = NULL) +
  facet_wrap(~panel) + coord_fixed() + labs(x = NULL, y = NULL, title = "Similarity between environments") +
  theme_paper() + theme(panel.grid = element_blank())
save_fig(f1, "fig01_environment_similarity.png", 9, 3.6)

# Fig 2: accuracy by validation design and model (points = environments)
f2 <- ggplot(by_env, aes(model, r, fill = model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) + geom_jitter(width = 0.1, size = 1.4, alpha = 0.7) +
  facet_wrap(~design, scales = "free_x") +
  scale_fill_manual(values = pal, guide = "none") +
  labs(x = NULL, y = "Predictive ability (r), one point per environment",
       title = "Prediction accuracy by validation design") +
  theme_paper() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_fig(f2, "fig02_accuracy_by_design.png", 10, 4.6)

# Fig 3: accuracy by environment for CV1 and CV2
f3 <- ggplot(subset(by_env, design != design_lab[["CV0_new_environment"]]),
             aes(env, r, colour = model, group = model)) +
  geom_line(linewidth = 0.6) + geom_point(size = 2) +
  facet_wrap(~design) + scale_colour_manual(values = pal, name = NULL) +
  labs(x = "Environment", y = "Predictive ability (r)", title = "Accuracy in each environment") + theme_paper()
save_fig(f3, "fig03_accuracy_by_environment.png", 9, 4)

# Fig 4: estimated vs realised genetic correlations
off <- upper.tri(cors$R_true)
f4d <- data.frame(truth = cors$R_true[off], estimate = cors$R_est[off])
f4 <- ggplot(f4d, aes(truth, estimate)) + geom_abline(linetype = 2) +
  geom_hline(yintercept = c(-1, 1), colour = "grey70") + geom_point(size = 2.2, colour = okabe_ito[1]) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "Realised genetic correlation (simulation truth)", y = "M4 estimated genetic correlation",
       title = "Unstructured model: correlation recovery",
       subtitle = "Estimates above 1 are not valid correlations\n(small-sample instability)") + theme_paper()
save_fig(f4, "fig04_correlation_recovery.png", 5.5, 4.4)

# Fig 5: new-environment prediction vs similarity of the held-out environment to the others
R <- cors$R_true; sim <- setNames((rowSums(R) - 1) / (ncol(R) - 1), rownames(R))
f5d <- subset(by_env, design == design_lab[["CV0_new_environment"]]); f5d$similarity <- sim[f5d$env]
f5 <- ggplot(f5d, aes(similarity, r, colour = model)) + geom_point(size = 2.5) +
  geom_text(aes(label = env), vjust = -1, size = 3, show.legend = FALSE) +
  scale_colour_manual(values = pal, name = NULL) +
  labs(x = "Mean genetic correlation with the other environments\n(simulation truth)",
       y = "Predictive ability (r)", title = "Predicting an unseen environment (CV0)") + theme_paper()
save_fig(f5, "fig05_new_environment.png", 6.5, 4.4)
# Fig 6: paired accuracy difference relative to the main-effect model M2
pd <- read.csv("results/tables/cv_paired_vs_M2.csv")
pd$design <- factor(design_lab[pd$design], levels = design_lab)
pd$model <- factor(model_labels[pd$model], levels = model_labels)
f6 <- ggplot(pd, aes(model, mean_diff_vs_M2, colour = model)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_pointrange(aes(ymin = mean_diff_vs_M2 - 2 * se, ymax = mean_diff_vs_M2 + 2 * se), size = 0.6) +
  facet_wrap(~design, scales = "free_x") + scale_colour_manual(values = pal, guide = "none") +
  labs(x = NULL, y = "Accuracy gain over M2 (r, mean \u00b1 2 SE)",
       title = "Paired comparison with the main-effect model",
       subtitle = "SEs from repeated CV are optimistic (folds share data); treat as a guide") +
  theme_paper() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_fig(f6, "fig06_gain_over_main_effect.png", 10, 4.6)
message("Saved 6 figures to results/figures/")
