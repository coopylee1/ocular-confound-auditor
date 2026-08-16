# =============================================================================
# run_revision2.R  --  E7b: does the gaze-trained subspace beat a rank-matched
# random projection, and if so WHERE?
#
# E7 (run_revision.R) found that at an 8-degree saccade amplitude every
# projection -- including a random rank-4 subspace -- recovers OPI ~ 0.43. That
# is not a flattering result for the instrument, and it needs explaining rather
# than hiding: at that amplitude the ocular signal (16 uV/deg * 8 deg) dwarfs
# both the neural source (1.8 uV) and the noise (3 uV), so almost any subspace
# retains it. The discriminating regime must therefore be the LOW-amplitude one
# -- which is also the regime the paper actually cares about, since sub-degree
# residual movements are what survive eye-movement rejection.
#
# E7b sweeps saccade amplitude x method, on the positive control (is the ocular
# signal still detected?) and on the covert-attention neural control (does the
# method stay null?). A useful instrument must do BOTH at the same amplitude.
# =============================================================================
suppressMessages({ library(parallel) })
source("R/auditor.R"); source("R/simulate.R")

OUT_T  <- "results/tables"
NPERM  <- 200
NSTUDY <- 12
NCORES <- max(1, min(detectCores() - 1, 10))
LAMBDA <- 1.0

src_full <- simulate_dataset("neutral", n_subjects = 20, trials_per = 200,
                             sacc_deg = 8, montage = "full", seed = 100)
AUD_full <- train_auditor(src_full$X, src_full$gaze, lambda = LAMBDA, subject = src_full$subj)
K_GAZE   <- AUD_full$rank

MONT    <- make_montage()
FRONTAL <- c("Fp1","Fpz","Fp2","AF7","AF3","AFz","AF4","AF8",
             "F7","F5","F3","F1","Fz","F2","F4","F6","F8")

proj_from_basis <- function(B) {
  q <- qr(B); U <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
  list(P = U %*% t(U), rank = q$rank)
}
proj_random <- function(p, k, seed) { set.seed(seed); proj_from_basis(matrix(rnorm(p * k), p, k)) }
proj_eog <- function(chans) {
  p <- length(chans); idx <- function(nm) match(nm, chans)
  v_h <- numeric(p); v_h[idx("F8")]  <-  1;  v_h[idx("F7")]  <- -1
  v_v <- numeric(p); v_v[idx("Fp1")] <- 0.5; v_v[idx("Fp2")] <- 0.5
  proj_from_basis(cbind(v_h, v_v))
}
proj_frontal <- function(chans) {
  keep <- intersect(FRONTAL, chans)
  B <- matrix(0, length(chans), length(keep))
  for (j in seq_along(keep)) B[match(keep[j], chans), j] <- 1
  proj_from_basis(B)
}

ci95 <- function(v) { v <- v[is.finite(v)]; qt(0.975, length(v) - 1) * sd(v) / sqrt(length(v)) }

cell <- function(seed, scenario, sacc_deg) {
  d <- simulate_dataset(scenario, n_subjects = 12, trials_per = 120,
                        coupling = 0.85, sacc_deg = sacc_deg,
                        montage = "full", seed = seed * 7919 + 13)
  chans <- d$channels; p <- length(chans)
  arms <- list(gaze    = list(P = AUD_full$P, rank = AUD_full$rank),
               random  = proj_random(p, K_GAZE, seed = seed * 31 + 7),
               eog     = proj_eog(chans),
               frontal = proj_frontal(chans))
  do.call(rbind, lapply(names(arms), function(nm) {
    o <- compute_opi(d$X, d$y, d$subj, arms[[nm]]$P, lambda = LAMBDA,
                     n_perm = NPERM, seed = seed)
    data.frame(seed = seed, scenario = scenario, sacc_deg = sacc_deg, method = nm,
               rank = arms[[nm]]$rank, opi = o$opi, auc_proj = o$auc_proj,
               auc_ceiling = o$auc_ceiling, p_value = o$p_value, row.names = NULL)
  }))
}

amps <- c(0.25, 0.5, 1, 2, 8)
cat("=== E7b: method x saccade amplitude ===\n")
E7b <- do.call(rbind, lapply(c("ocular_corr", "covert"), function(sc)
  do.call(rbind, lapply(amps, function(a) {
    cat(sprintf("  %s  amp=%.2f deg\n", sc, a))
    do.call(rbind, mclapply(seq_len(NSTUDY), function(s) cell(s, sc, a), mc.cores = NCORES))
  }))))

E7bs <- do.call(rbind, lapply(split(E7b, list(E7b$scenario, E7b$sacc_deg, E7b$method), drop = TRUE),
  function(g) data.frame(scenario = g$scenario[1], sacc_deg = g$sacc_deg[1],
                         method = g$method[1], rank = g$rank[1],
                         opi_m = mean(g$opi), opi_sd = sd(g$opi), opi_ci = ci95(g$opi),
                         proj_m = mean(g$auc_proj), ceil_m = mean(g$auc_ceiling),
                         frac_sig = mean(g$p_value < 0.05), row.names = NULL)))
E7bs <- E7bs[order(E7bs$scenario, E7bs$sacc_deg, E7bs$method), ]
print(E7bs[, c("scenario","sacc_deg","method","opi_m","opi_ci","frac_sig")],
      digits = 3, row.names = FALSE)

cat("\n--- Wide view: OPI by amplitude and method ---\n")
for (sc in unique(E7bs$scenario)) {
  cat(sprintf("\n%s:\n", sc))
  sub <- E7bs[E7bs$scenario == sc, ]
  w <- reshape(sub[, c("sacc_deg","method","opi_m")], idvar = "sacc_deg",
               timevar = "method", direction = "wide")
  print(w, digits = 3, row.names = FALSE)
}

write.csv(E7b,  file.path(OUT_T, "E7b_amplitude_methods_raw.csv"),     row.names = FALSE)
write.csv(E7bs, file.path(OUT_T, "E7b_amplitude_methods_summary.csv"), row.names = FALSE)
cat("\n=== run_revision2.R DONE ===\n")
