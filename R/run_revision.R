# =============================================================================
# run_revision.R  --  Analyses added for the IEEE Healthcom 2026 camera-ready,
# answering the reviewers' specific technical requests.
#
#   E7  Baselines (R1): the OPI computed under FOUR projections of stated rank --
#       (a) the proposed gaze-trained subspace, (b) a RANDOM subspace matched in
#       rank, (c) an EOG-proxy built from channel arithmetic only, and (d) a
#       frontal-channel audit. The discriminating question is not which is most
#       sensitive on positive controls, but which stays NULL on the neural
#       negative controls.
#
#   E8  Per-output reliability (R1): reconstruction reported SEPARATELY for
#       horizontal gaze, vertical gaze and the saccade channel (R^2, continuous)
#       and for blinks (AUC, because the blink target is binary) -- replacing the
#       single unexplained aggregate.
#
#   E9  Haufe activation pattern (R1): the neural-overlap score recomputed with
#       a = Cov(X) w, a valid activation pattern, instead of the raw
#       discriminative weight vector w.
#
#   E10 Held-out threshold calibration (R1): thresholds fit on the full-64
#       montage grid ONLY, then applied to the held-out ERP-CORE-30 grid.
#
# Every number is produced by running this file; nothing is asserted.
# =============================================================================
suppressMessages({ library(parallel); library(jsonlite) })
source("R/auditor.R"); source("R/simulate.R")

OUT_T  <- "results/tables"
dir.create(OUT_T, showWarnings = FALSE, recursive = TRUE)
NPERM  <- 200
NSTUDY <- 12
NCORES <- max(1, min(detectCores() - 1, 10))
LAMBDA <- 1.0

set.seed(20260815)

cat("=== Frozen auditors (same source as run_simulation.R, seed 100) ===\n")
src_full <- simulate_dataset("neutral", n_subjects = 20, trials_per = 200,
                             sacc_deg = 8, montage = "full",    seed = 100)
src_erp  <- simulate_dataset("neutral", n_subjects = 20, trials_per = 200,
                             sacc_deg = 8, montage = "erpcore", seed = 100)
AUD <- list(
  full    = train_auditor(src_full$X, src_full$gaze, lambda = LAMBDA, subject = src_full$subj),
  erpcore = train_auditor(src_erp$X,  src_erp$gaze,  lambda = LAMBDA)
)
K_GAZE <- AUD$full$rank
cat(sprintf("  gaze-predictive subspace rank: full(64ch)=%d  erpcore(30ch)=%d\n\n",
            AUD$full$rank, AUD$erpcore$rank))

MONT     <- make_montage()
CH_FULL  <- MONT$full
FRONTAL  <- c("Fp1","Fpz","Fp2","AF7","AF3","AFz","AF4","AF8",
              "F7","F5","F3","F1","Fz","F2","F4","F6","F8")

# ---------------------------------------------------------------------------
# Comparator projections. Each is an orthogonal projector on the SAME channel
# space, so the only thing that differs between arms is WHICH subspace is used.
# ---------------------------------------------------------------------------
proj_from_basis <- function(B) {                      # B: p x k, columns span
  q <- qr(B); U <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
  list(P = U %*% t(U), rank = q$rank)
}

# (b) random subspace, matched in rank to the gaze-trained operator
proj_random <- function(p, k, seed) {
  set.seed(seed)
  proj_from_basis(matrix(rnorm(p * k), p, k))
}

# (c) EOG proxy: bipolar channel arithmetic only -- no gaze training whatsoever.
#     HEOG ~ F8 - F7 (horizontal, antisymmetric), VEOG ~ (Fp1 + Fp2)/2 (vertical).
proj_eog <- function(chans) {
  p <- length(chans); idx <- function(nm) match(nm, chans)
  v_h <- numeric(p); v_h[idx("F8")]  <-  1; v_h[idx("F7")]  <- -1
  v_v <- numeric(p); v_v[idx("Fp1")] <- 0.5; v_v[idx("Fp2")] <- 0.5
  proj_from_basis(cbind(v_h, v_v))
}

# (d) frontal-channel audit: coordinate subspace spanned by frontal electrodes
proj_frontal <- function(chans) {
  keep <- intersect(FRONTAL, chans)
  B <- matrix(0, length(chans), length(keep))
  for (j in seq_along(keep)) B[match(keep[j], chans), j] <- 1
  proj_from_basis(B)
}

# ---------------------------------------------------------------------------
# E7: OPI under each projection, across positive and negative controls.
# ---------------------------------------------------------------------------
cat("=== E7: rank-matched and simple-alternative baselines ===\n")

e7_cell <- function(seed, scenario, coupling = 0.85) {
  d <- simulate_dataset(scenario, n_subjects = 12, trials_per = 120,
                        coupling = coupling, sacc_deg = 8,
                        montage = "full", seed = seed * 7919 + 13)
  chans <- d$channels; p <- length(chans)
  arms <- list(
    gaze    = list(P = AUD$full$P,                  rank = AUD$full$rank),
    random  = proj_random(p, K_GAZE, seed = seed * 31 + 7),
    eog     = proj_eog(chans),
    frontal = proj_frontal(chans)
  )
  do.call(rbind, lapply(names(arms), function(nm) {
    o <- compute_opi(d$X, d$y, d$subj, arms[[nm]]$P, lambda = LAMBDA,
                     n_perm = NPERM, seed = seed)
    data.frame(seed = seed, scenario = scenario, method = nm,
               rank = arms[[nm]]$rank, opi = o$opi, auc_proj = o$auc_proj,
               auc_ceiling = o$auc_ceiling, p_value = o$p_value,
               row.names = NULL)
  }))
}

scen7 <- c("ocular_corr", "blink_corr", "covert", "motor")
E7 <- do.call(rbind, lapply(scen7, function(sc) {
  cat(sprintf("  scenario: %s\n", sc))
  rows <- mclapply(seq_len(NSTUDY), function(s) e7_cell(s, sc), mc.cores = NCORES)
  do.call(rbind, rows)
}))

ci95 <- function(v) { v <- v[is.finite(v)]; qt(0.975, length(v) - 1) * sd(v) / sqrt(length(v)) }
E7s <- do.call(rbind, lapply(split(E7, list(E7$scenario, E7$method), drop = TRUE), function(g)
  data.frame(scenario = g$scenario[1], method = g$method[1], rank = g$rank[1],
             opi_m = mean(g$opi), opi_sd = sd(g$opi), opi_ci = ci95(g$opi),
             proj_m = mean(g$auc_proj), ceil_m = mean(g$auc_ceiling),
             frac_sig = mean(g$p_value < 0.05), row.names = NULL)))
E7s <- E7s[order(E7s$scenario, E7s$method), ]
print(E7s, digits = 3, row.names = FALSE)
write.csv(E7,  file.path(OUT_T, "E7_baselines_raw.csv"),     row.names = FALSE)
write.csv(E7s, file.path(OUT_T, "E7_baselines_summary.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# E8: per-output reconstruction. R^2 for the three continuous ocular targets;
# AUC for the binary blink target (an R^2 on a 0/1 outcome is not interpretable).
# ---------------------------------------------------------------------------
cat("\n=== E8: per-output reliability (gaze_x, gaze_y, saccade: R2; blink: AUC) ===\n")

reliability_per_output <- function(X, Y, subj, lambda = LAMBDA) {
  Y <- as.matrix(Y); subs <- unique(subj)
  pred <- matrix(NA_real_, nrow(X), ncol(Y))
  for (s in subs) {
    tr <- subj != s; te <- subj == s
    fit <- ridge_fit(X[tr, , drop = FALSE], Y[tr, , drop = FALSE], lambda)
    pred[te, ] <- ridge_predict(fit, X[te, , drop = FALSE])
  }
  r2 <- gaze_r2(Y, pred)
  blink_col <- which(colnames(Y) == "blink")
  blink_auc <- if (length(blink_col) == 1 && length(unique(Y[, blink_col])) == 2)
    auc_score(pred[, blink_col], Y[, blink_col]) else NA_real_
  list(r2_gx = unname(r2[1]), r2_gy = unname(r2[2]), r2_sacc = unname(r2[3]),
       blink_auc = blink_auc, mean_r2_all = mean(r2, na.rm = TRUE),
       mean_r2_cont = mean(r2[1:3], na.rm = TRUE))
}

e8_cell <- function(seed, montage, sacc_deg) {
  d <- simulate_dataset("ocular_corr", n_subjects = 12, trials_per = 120,
                        coupling = 0.85, sacc_deg = sacc_deg,
                        montage = montage, seed = seed * 7919 + 13)
  r <- reliability_per_output(d$X, d$gaze, d$subj)
  data.frame(seed = seed, montage = montage, sacc_deg = sacc_deg,
             r2_gx = r$r2_gx, r2_gy = r$r2_gy, r2_sacc = r$r2_sacc,
             blink_auc = r$blink_auc, mean_r2_all = r$mean_r2_all,
             mean_r2_cont = r$mean_r2_cont, row.names = NULL)
}

amps <- c(0.25, 0.5, 1, 2, 4, 8, 12)
E8 <- do.call(rbind, lapply(amps, function(a) do.call(rbind, lapply(c("full","erpcore"), function(m) {
  rows <- mclapply(seq_len(NSTUDY), function(s) e8_cell(s, m, a), mc.cores = NCORES)
  do.call(rbind, rows)
}))))
E8s <- do.call(rbind, lapply(split(E8, list(E8$montage, E8$sacc_deg), drop = TRUE), function(g)
  data.frame(montage = g$montage[1], sacc_deg = g$sacc_deg[1],
             r2_gx = mean(g$r2_gx),   r2_gx_sd = sd(g$r2_gx),
             r2_gy = mean(g$r2_gy),   r2_gy_sd = sd(g$r2_gy),
             r2_sacc = mean(g$r2_sacc), r2_sacc_sd = sd(g$r2_sacc),
             blink_auc = mean(g$blink_auc), blink_auc_sd = sd(g$blink_auc),
             mean_r2_all = mean(g$mean_r2_all), mean_r2_cont = mean(g$mean_r2_cont),
             row.names = NULL)))
E8s <- E8s[order(E8s$montage, E8s$sacc_deg), ]
print(E8s[, c("montage","sacc_deg","r2_gx","r2_gy","r2_sacc","blink_auc","mean_r2_all")],
      digits = 3, row.names = FALSE)
write.csv(E8,  file.path(OUT_T, "E8_per_output_raw.csv"),     row.names = FALSE)
write.csv(E8s, file.path(OUT_T, "E8_per_output_summary.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# E9: overlap computed on the Haufe activation pattern a = Cov(X) w rather than
# on the discriminative weight vector w. R1 is right that w is not an activation
# pattern; we report both so the change is auditable.
# ---------------------------------------------------------------------------
cat("\n=== E9: neural leakage vs overlap, weights vs Haufe activation pattern ===\n")

C_FULL   <- MONT$coords[CH_FULL, , drop = FALSE]
topo_all <- make_topographies(list(coords = C_FULL, eye_left = MONT$eye_left,
                                   eye_right = MONT$eye_right))
P_FULL   <- AUD$full$P
overlap_of <- function(t) sqrt(sum((P_FULL %*% t)^2)) / sqrt(sum(t^2))

neural_contrast_topo <- function(ypos)
  .prox(C_FULL, c(0.50, ypos), 0.28) - .prox(C_FULL, c(-0.50, ypos), 0.28)

gen_neural_only <- function(ypos, seed, n_subjects = 10, trials_per = 120, neural_uV = 2.2) {
  set.seed(seed)
  N <- n_subjects * trials_per; subj <- rep(seq_len(n_subjects), each = trials_per)
  y <- rbinom(N, 1, 0.5); side <- ifelse(y == 1, 1, -1)
  X <- .bg_noise(N, nrow(C_FULL), C_FULL, sd_smooth = 3, sd_white = 3)
  X <- X + outer(neural_uV * side, neural_contrast_topo(ypos))
  gx <- rnorm(N, 0, 8) * sample(c(-1, 1), N, TRUE)
  X <- X + outer(gx * 16, topo_all$heog) + outer(rnorm(N, 0, 2.5) * 6, topo_all$veog)
  list(X = X, y = y, subj = subj)
}

# Empirical contrast direction as a decoder would estimate it, and its Haufe
# transform: a = Cov(X) w  (Haufe et al. 2014). Both are unit-normalised before
# the overlap is taken, so s is scale-free.
fit_w_and_haufe <- function(X, y) {
  fit <- ridge_fit(X, matrix(2 * y - 1, ncol = 1), lambda = LAMBDA)
  w   <- as.numeric(fit$W)
  Xc  <- sweep(X, 2L, colMeans(X), "-")
  a   <- as.numeric(crossprod(Xc) %*% w) / (nrow(X) - 1)
  list(w = w, a = a)
}

ypositions <- seq(0.55, -0.90, length.out = 10)
E9 <- do.call(rbind, mclapply(ypositions, function(yp) {
  s_true <- overlap_of(neural_contrast_topo(yp))     # overlap of the TRUE source topography
  per <- lapply(1:8, function(s) {
    d  <- gen_neural_only(yp, seed = s * 101 + 5)
    wa <- fit_w_and_haufe(d$X, d$y)
    o  <- compute_opi(d$X, d$y, d$subj, P_FULL, lambda = LAMBDA, n_perm = 150, seed = s)
    c(opi = o$opi, s_w = overlap_of(wa$w), s_a = overlap_of(wa$a))
  })
  M <- do.call(rbind, per)
  data.frame(ypos = yp, s_true = s_true,
             s_weights = mean(M[, "s_w"]), s_weights_sd = sd(M[, "s_w"]),
             s_haufe   = mean(M[, "s_a"]), s_haufe_sd   = sd(M[, "s_a"]),
             opi_leak  = mean(M[, "opi"]), opi_sd = sd(M[, "opi"]),
             opi_ci    = ci95(M[, "opi"]), row.names = NULL)
}, mc.cores = max(1, min(detectCores() - 1, 8))))
E9 <- E9[order(E9$s_haufe), ]
print(E9[, c("s_true","s_weights","s_haufe","opi_leak","opi_sd")], digits = 3, row.names = FALSE)
rho_w <- cor(E9$s_weights, E9$opi_leak, method = "spearman")
rho_a <- cor(E9$s_haufe,   E9$opi_leak, method = "spearman")
rho_t <- cor(E9$s_true,    E9$opi_leak, method = "spearman")
cat(sprintf("\nSpearman rho(overlap, leakage):  weights=%.3f   Haufe=%.3f   true-topo=%.3f\n",
            rho_w, rho_a, rho_t))
write.csv(E9, file.path(OUT_T, "E9_haufe_overlap.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# E10: threshold calibration on a HELD-OUT grid. Fit thresholds using the
# full-64 montage cells only; then apply them, unchanged, to the ERP-CORE-30
# cells that played no part in choosing them.
# ---------------------------------------------------------------------------
cat("\n=== E10: held-out reliability-threshold calibration ===\n")
E3s_path <- file.path(OUT_T, "E3_amplitude_summary.csv")
if (file.exists(E3s_path)) {
  E3s <- read.csv(E3s_path)
  fit_grid  <- E3s[E3s$montage == "full",    ]
  test_grid <- E3s[E3s$montage == "erpcore", ]
  pick <- function(g, tol) { ok <- g$r2_m[abs(g$bias_m) <= tol]; if (length(ok)) min(ok) else NA_real_ }
  t_conf_ho <- pick(fit_grid, 0.05)
  t_expl_ho <- pick(fit_grid, 0.10)
  cat(sprintf("  fit on full-64 grid only:  t_conf=%.3f  t_expl=%.3f  (rounded: %.2f / %.2f)\n",
              t_conf_ho, t_expl_ho, round(t_conf_ho, 2), round(t_expl_ho, 2)))
  test_grid$tier_heldout <- mapply(reliability_tier, test_grid$r2_m,
                                   MoreArgs = list(t_conf = round(t_conf_ho, 2),
                                                   t_expl = round(t_expl_ho, 2)))
  print(test_grid[, c("montage","sacc_deg","r2_m","bias_m","tier_heldout")],
        digits = 3, row.names = FALSE)
  # Does the held-out gate actually separate acceptable from unacceptable bias?
  ok_bias <- abs(test_grid$bias_m) <= 0.05
  cat(sprintf("  held-out cells tiered confirmatory: %d/%d; of those, |bias|<=0.05 in %d\n",
              sum(test_grid$tier_heldout == "confirmatory"), nrow(test_grid),
              sum(test_grid$tier_heldout == "confirmatory" & ok_bias)))
  write_json(list(t_conf_heldout = t_conf_ho, t_expl_heldout = t_expl_ho,
                  t_conf_reported = round(t_conf_ho, 2), t_expl_reported = round(t_expl_ho, 2),
                  fit_montage = "full64", test_montage = "erpcore30",
                  test_grid = test_grid),
             file.path(OUT_T, "E10_threshold_heldout.json"),
             digits = 5, auto_unbox = TRUE, pretty = TRUE)
  write.csv(test_grid, file.path(OUT_T, "E10_threshold_heldout.csv"), row.names = FALSE)
} else {
  cat("  SKIPPED: results/tables/E3_amplitude_summary.csv not found; run run_simulation.R first.\n")
}

cat("\n=== run_revision.R DONE ===\n")
saveRDS(list(E7 = E7, E7s = E7s, E8 = E8, E8s = E8s, E9 = E9,
             rho = c(weights = rho_w, haufe = rho_a, true = rho_t)),
        file.path(OUT_T, "revision_results.rds"))
