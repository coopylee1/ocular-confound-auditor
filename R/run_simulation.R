# =============================================================================
# run_simulation.R  --  the in-silico validation of the Ocular-Confound Auditor
#
# Produces (saved to results/tables/ and results/figures/):
#   E1  H3 sensitivity battery     - OPI high on condition-correlated ocular;
#                                    ~null on covert-attention, motor, neutral,
#                                    label-shuffle. (positive + negative controls)
#   E2  H3 graded coupling sweep   - OPI rises monotonically with gaze<->label
#                                    coupling.
#   E3  Method 5 amplitude / bias  - across the saccade-amplitude range
#       envelope                     (fixational -> large), measure auditor OPI
#                                    vs ground-truth-gaze OPI (signed two-sided
#                                    bias) and gaze-reconstruction reliability,
#                                    on dense vs ERP-CORE-30 montage.
#   E4  Method 3 reliability gate  - empirical bias-vs-fidelity curve -> the
#                                    confirmatory/exploratory/inconclusive
#                                    thresholds; demonstrate "inconclusive" for
#                                    a fixational signal on the sparse montage.
#   E5  H2-analog headline         - simulated ERP-CORE-like N2pc audit with a
#                                    small surviving sub-0.2-deg ocular residual.
#
# Every number is produced by actually running the code below; nothing is
# asserted. RNG seeds are fixed for exact reproducibility.
# =============================================================================
suppressMessages({ library(parallel); library(jsonlite) })
source("R/auditor.R"); source("R/simulate.R")

OUT_T <- "results/tables"; OUT_F <- "results/figures"
dir.create(OUT_T, showWarnings = FALSE, recursive = TRUE)
NPERM   <- 200
NSTUDY  <- 12                              # independent simulated "studies" per cell
NCORES  <- max(1, min(detectCores() - 1, 10))
LAMBDA  <- 1.0

cat("=== Training the frozen auditor on a large-saccade source (EEGEyeNet-like) ===\n")
# EEGEyeNet has ONLY large cued saccades (fact-checked) -> train at sacc_deg=8.
src_full <- simulate_dataset("neutral", n_subjects = 20, trials_per = 200,
                              sacc_deg = 8, montage = "full",   seed = 100)
src_erp   <- simulate_dataset("neutral", n_subjects = 20, trials_per = 200,
                              sacc_deg = 8, montage = "erpcore", seed = 100)
AUD <- list(
  full    = train_auditor(src_full$X, src_full$gaze, lambda = LAMBDA, subject = src_full$subj),
  erpcore = train_auditor(src_erp$X,   src_erp$gaze,   lambda = LAMBDA)
)
cat(sprintf("  auditor subspace rank: full(64ch)=%d  erpcore(30ch)=%d\n\n",
            AUD$full$rank, AUD$erpcore$rank))

# Helper: one simulated study -> a result row.
run_cell <- function(seed, scenario, coupling = 0.8, sacc_deg = 8,
                     montage = "full", shuffle_labels = FALSE,
                     n_subjects = 12, trials_per = 120, blink_rate = 0.12) {
  d <- simulate_dataset(scenario, n_subjects = n_subjects, trials_per = trials_per,
                        coupling = coupling, sacc_deg = sacc_deg,
                        blink_rate = blink_rate, montage = montage,
                        seed = seed * 7919 + 13)
  y <- d$y
  if (shuffle_labels) {                         # destroy label assoc within subject
    for (s in unique(d$subj)) { i <- which(d$subj == s); y[i] <- sample(y[i]) }
  }
  aud <- AUD[[montage]]
  rel <- auditor_reliability(d$X, d$gaze, d$subj, lambda = LAMBDA)
  o   <- compute_opi(d$X, y, d$subj, aud$P, lambda = LAMBDA,
                     n_perm = NPERM, seed = seed)
  ot  <- compute_opi_from_truth(d$gaze, y, d$subj, lambda = LAMBDA,
                                n_perm = NPERM, seed = seed)
  data.frame(seed = seed, scenario = scenario, montage = montage,
             coupling = coupling, sacc_deg = sacc_deg,
             shuffle = shuffle_labels,
             opi = o$opi, auc_proj = o$auc_proj, auc_ceiling = o$auc_ceiling,
             null_mean = o$null_mean, p_value = o$p_value,
             opi_truth = ot$opi_truth, bias = o$opi - ot$opi_truth,
             gaze_r2 = rel$mean_r2,
             r2_gx = rel$per_output_r2[1], r2_blink = rel$per_output_r2[4])
}

mc <- function(seeds, ...) {
  args <- list(...)
  rows <- mclapply(seeds, function(s) do.call(run_cell, c(list(seed = s), args)),
                   mc.cores = NCORES)
  do.call(rbind, rows)
}
summ <- function(df, by) {
  agg <- function(v) c(mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE))
  s <- do.call(rbind, lapply(split(df, df[by]), function(g) {
    data.frame(g[1, by, drop = FALSE],
               opi_m = mean(g$opi), opi_sd = sd(g$opi),
               proj_m = mean(g$auc_proj), ceil_m = mean(g$auc_ceiling),
               bias_m = mean(g$bias), bias_sd = sd(g$bias),
               r2_m = mean(g$gaze_r2),
               p_med = median(g$p_value),
               frac_sig = mean(g$p_value < 0.05), row.names = NULL)
  }))
  s
}

# ---------------------------------------------------------------- E1: H3 battery
cat("=== E1: H3 sensitivity battery (positive + negative controls) ===\n")
seeds <- seq_len(NSTUDY)
E1 <- rbind(
  mc(seeds, scenario = "ocular_corr", coupling = 0.85, sacc_deg = 8, montage = "full"),
  mc(seeds, scenario = "blink_corr",  coupling = 0.85, sacc_deg = 8, montage = "full"),
  mc(seeds, scenario = "covert",                       sacc_deg = 8, montage = "full"),
  mc(seeds, scenario = "motor",                        sacc_deg = 8, montage = "full"),
  mc(seeds, scenario = "neutral",                      sacc_deg = 8, montage = "full"),
  mc(seeds, scenario = "ocular_corr", coupling = 0.85, sacc_deg = 8, montage = "full",
     shuffle_labels = TRUE)
)
E1$scenario[E1$shuffle] <- "label_shuffle"
E1s <- summ(E1, "scenario")
print(E1s[order(-E1s$opi_m), c("scenario","opi_m","opi_sd","proj_m","ceil_m","p_med","frac_sig")],
      digits = 3, row.names = FALSE)
write.csv(E1,  file.path(OUT_T, "E1_battery_raw.csv"), row.names = FALSE)
write.csv(E1s, file.path(OUT_T, "E1_battery_summary.csv"), row.names = FALSE)

# ----------------------------------------------------- E2: graded coupling sweep
cat("\n=== E2: H3 graded coupling sweep (OPI vs gaze-label coupling) ===\n")
E2 <- do.call(rbind, lapply(c(0,0.2,0.4,0.6,0.8,1.0), function(rho)
  mc(seeds, scenario = "ocular_corr", coupling = rho, sacc_deg = 8, montage = "full")))
E2s <- summ(E2, "coupling")
print(E2s[, c("coupling","opi_m","opi_sd","proj_m","frac_sig")], digits = 3, row.names = FALSE)
write.csv(E2,  file.path(OUT_T, "E2_coupling_raw.csv"), row.names = FALSE)
write.csv(E2s, file.path(OUT_T, "E2_coupling_summary.csv"), row.names = FALSE)

# -------------------------------------- E3/E4: amplitude, bias envelope, gate
cat("\n=== E3/E4: amplitude sweep, bias envelope & reliability (full-64 vs erpcore-30) ===\n")
amps <- c(0.25, 0.5, 1, 2, 4, 8, 12)
E3 <- do.call(rbind, lapply(amps, function(a) rbind(
  mc(seeds, scenario = "ocular_corr", coupling = 0.85, sacc_deg = a, montage = "full"),
  mc(seeds, scenario = "ocular_corr", coupling = 0.85, sacc_deg = a, montage = "erpcore")
)))
E3s <- summ(E3, c("montage","sacc_deg"))
E3s <- E3s[order(E3s$montage, E3s$sacc_deg), ]
print(E3s[, c("montage","sacc_deg","r2_m","opi_m","bias_m","bias_sd")], digits = 3, row.names = FALSE)
write.csv(E3,  file.path(OUT_T, "E3_amplitude_raw.csv"), row.names = FALSE)
write.csv(E3s, file.path(OUT_T, "E3_amplitude_summary.csv"), row.names = FALSE)

# Empirical reliability thresholds (Method 3): bias-vs-fidelity curve.
# Choose t_conf = smallest gaze-R^2 with mean |bias| <= 0.05 (confirmatory),
#        t_expl = smallest gaze-R^2 with mean |bias| <= 0.10 (exploratory).
fid <- E3s[order(E3s$r2_m), ]
abias <- abs(fid$bias_m)
t_conf <- { ok <- fid$r2_m[abias <= 0.05]; if (length(ok)) min(ok) else NA }
t_expl <- { ok <- fid$r2_m[abias <= 0.10]; if (length(ok)) min(ok) else NA }
cat(sprintf("\nCalibrated reliability thresholds:  t_conf(R2)=%.3f  t_expl(R2)=%.3f\n",
            t_conf, t_expl))
calib <- list(t_conf = t_conf, t_expl = t_expl,
              bias_fidelity = data.frame(r2 = fid$r2_m, abs_bias = abias,
                                         montage = fid$montage, sacc_deg = fid$sacc_deg))
write_json(calib, file.path(OUT_T, "calibration.json"), digits = 5, auto_unbox = TRUE, pretty = TRUE)

# Tier each amplitude/montage cell with the calibrated gate.
E3s$tier <- mapply(reliability_tier, E3s$r2_m,
                   MoreArgs = list(t_conf = t_conf, t_expl = t_expl))
write.csv(E3s, file.path(OUT_T, "E3_amplitude_summary.csv"), row.names = FALSE)
cat("\nReliability tiers across the amplitude x montage grid:\n")
print(E3s[, c("montage","sacc_deg","r2_m","opi_m","tier")], digits = 3, row.names = FALSE)

# ---------------------------------------------------- E5: H2-analog headline
cat("\n=== E5: H2-analog -- simulated ERP-CORE N2pc audit with sub-0.2-deg residual ===\n")
# ERP CORE controls eye movements: N2pc trials with >0.2 deg HEOG are rejected.
# We emulate the SURVIVING residual: a tiny condition-correlated fixational
# movement (~0.15 deg) on top of the genuine posterior N2pc, sparse montage.
E5 <- do.call(rbind, lapply(seeds, function(s) {
  d  <- simulate_dataset("covert", n_subjects = 12, trials_per = 160,
                         sacc_deg = 0.15, montage = "erpcore", seed = s * 7919 + 13)
  # inject a small condition-correlated horizontal residual (<0.2 deg)
  resid <- 0.15 * d$side * rbinom(nrow(d$X), 1, 0.7)
  mont  <- make_montage(); topo <- make_topographies(list(coords = mont$coords[d$channels,],
            eye_left = mont$eye_left, eye_right = mont$eye_right))
  d$X   <- d$X + outer(resid * 16, topo$heog)
  d$gaze[,"gx"] <- d$gaze[,"gx"] + resid
  aud <- AUD$erpcore
  rel <- auditor_reliability(d$X, d$gaze, d$subj, lambda = LAMBDA)
  o   <- compute_opi(d$X, d$y, d$subj, aud$P, lambda = LAMBDA, n_perm = NPERM, seed = s)
  data.frame(seed = s, opi = o$opi, auc_proj = o$auc_proj, auc_ceiling = o$auc_ceiling,
             p_value = o$p_value, gaze_r2 = rel$mean_r2,
             tier = reliability_tier(rel$mean_r2, t_conf, t_expl))
}))
cat(sprintf("  mean OPI=%+.3f (sd %.3f)  ceilingAUC=%.3f  median p=%.3f  frac p<.05=%.2f\n",
            mean(E5$opi), sd(E5$opi), mean(E5$auc_ceiling), median(E5$p_value),
            mean(E5$p_value < 0.05)))
cat(sprintf("  reliability tier (modal): %s   (mean gaze R2=%.3f)\n",
            names(sort(table(E5$tier), decreasing = TRUE))[1], mean(E5$gaze_r2)))
write.csv(E5, file.path(OUT_T, "E5_h2analog.csv"), row.names = FALSE)

cat("\n=== DONE. Tables written to results/tables/. ===\n")
saveRDS(list(E1=E1,E2=E2,E3=E3,E5=E5,calib=calib,AUD_rank=c(AUD$full$rank,AUD$erpcore$rank)),
        file.path(OUT_T, "sim_results.rds"))
