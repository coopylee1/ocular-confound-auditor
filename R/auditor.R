# =============================================================================
# auditor.R  --  The Ocular-Confound Auditor (core library)
#
# Implements, faithfully, the LOCKED definitions from the paper proposal:
#   * Method 1: the ocular-subspace projection operator is the ORTHOGONAL
#               projection onto the span of the frozen LINEAR gaze/ocular-state
#               readout weights.  (A linear ridge readout is exactly this
#               object -- no deep network is needed to instantiate the locked
#               operator.)
#   * Method 4: the Ocular Predictability Index (OPI) is the AUC of a
#               paradigm-matched cognitive decoder applied to the locked
#               ocular-only projection, minus a within-subject label-permutation
#               null, with full-EEG AUC as the ceiling.
#   * Method 3: an empirically-derived Auditor Reliability Score gates the
#               result into confirmatory / exploratory / inconclusive tiers.
#   * Method 5: a MEASURED two-sided bias envelope -- auditor OPI vs the OPI
#               computed from ground-truth gaze -- (see run_simulation.R).
#
# Pure functions only; no file I/O here. Everything is base R + the rank-based
# AUC implemented below (cross-checked against pROC in the test harness).
# =============================================================================

# ---- Area under the ROC curve (Mann-Whitney U form, ties handled by rank) ----
auc_score <- function(scores, labels) {
  labels <- as.integer(labels)
  pos <- scores[labels == 1L]
  neg <- scores[labels == 0L]
  n1 <- length(pos); n0 <- length(neg)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(c(pos, neg))                       # midranks -> tie-correct
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# ---- Ridge regression (multi-output), scale-invariant penalty -----------------
# W minimises ||Yc - Xc W||^2 + lambda * (tr(Xc'Xc)/p) * ||W||^2
# X: n x p, Y: n x q.  Returns weights plus the centring constants.
ridge_fit <- function(X, Y, lambda = 1.0) {
  X <- as.matrix(X); Y <- as.matrix(Y)
  mx <- colMeans(X); my <- colMeans(Y)
  Xc <- sweep(X, 2L, mx, "-")
  Yc <- sweep(Y, 2L, my, "-")
  p  <- ncol(Xc)
  XtX <- crossprod(Xc)
  ridge <- lambda * (sum(diag(XtX)) / p)       # dimensionless lambda
  W <- solve(XtX + ridge * diag(p), crossprod(Xc, Yc))
  list(W = W, mx = mx, my = my)
}

ridge_predict <- function(fit, Xnew) {
  Xc <- sweep(as.matrix(Xnew), 2L, fit$mx, "-")
  sweep(Xc %*% fit$W, 2L, fit$my, "+")
}

# =============================================================================
# THE LOCKED OPERATOR (Method 1)
# Orthogonal projector onto the column space (span) of the readout weights W.
# W is p x q (channels/features x ocular outputs). The ocular subspace is
# span(columns of W); its orthogonal projector in feature space is U U^T, where
# U holds the left singular vectors of W spanning that column space.
# =============================================================================
ocular_projection <- function(W, tol = 1e-8) {
  W <- as.matrix(W)
  s <- svd(W)
  k <- if (length(s$d) == 0L) 0L else sum(s$d > tol * s$d[1])
  if (k == 0L) return(matrix(0, nrow(W), nrow(W)))
  U <- s$u[, seq_len(k), drop = FALSE]
  list(P = U %*% t(U), rank = k, basis = U)
}

# =============================================================================
# Train the frozen auditor on SOURCE data (EEG features X_src + gaze/ocular
# ground-truth targets Y_src). Returns the readout, the locked projector P, and
# the subspace rank. This is the object that gets "frozen and shipped".
# =============================================================================
train_auditor <- function(X_src, Y_src, lambda = 1.0, tol = 1e-8, subject = NULL) {
  fit <- ridge_fit(X_src, Y_src, lambda = lambda)
  proj <- ocular_projection(fit$W, tol = tol)
  # Reliability is a property of the (montage-restricted) auditor measured on the
  # SOURCE -- it needs no target eye-tracker, matching the no-eye-tracker premise.
  rel <- if (!is.null(subject)) auditor_reliability(X_src, Y_src, subject, lambda)$mean_r2 else NA_real_
  structure(list(readout = fit, P = proj$P, rank = proj$rank,
                 basis = proj$basis, n_outputs = ncol(as.matrix(Y_src)),
                 reliability_r2 = rel, channels = colnames(as.matrix(X_src))),
            class = "ocular_auditor")
}

# =============================================================================
# THE SHIPPED TOOL (Method 7)
# audit_dataset(): EEG matrix + labels + subject + a frozen auditor  ->  the
# Ocular Predictability Index with its two-sided envelope and trust tier.
# CPU-runnable; no eye-tracker required on the target.
# =============================================================================
audit_dataset <- function(X, labels, subject, auditor,
                          t_conf = 0.50, t_expl = 0.20, bias_margin = 0.05,
                          lambda = 1.0, n_perm = 200, seed = 1) {
  stopifnot(inherits(auditor, "ocular_auditor"))
  o <- compute_opi(X, labels, subject, auditor$P, lambda = lambda,
                   n_perm = n_perm, seed = seed)
  tier <- reliability_tier(auditor$reliability_r2, t_conf, t_expl)
  out <- list(
    opi          = o$opi,                                   # the index
    opi_envelope = c(lower = o$opi - bias_margin,           # two-sided bias envelope
                     upper = o$opi + bias_margin),
    ceiling_auc  = o$auc_ceiling,
    projected_auc= o$auc_proj,
    null_mean    = o$null_mean,
    p_value      = o$p_value,
    reliability_r2 = auditor$reliability_r2,
    trust_tier   = tier,                                    # confirmatory/exploratory/inconclusive
    verdict      = if (tier == "inconclusive")
                     "INCONCLUSIVE: auditor cannot reconstruct gaze in this montage/regime"
                   else if (o$opi - bias_margin > 0)
                     sprintf("ocular leakage detected (OPI %.3f > margin %.3f), tier=%s",
                             o$opi, bias_margin, tier)
                   else
                     sprintf("no ocular leakage above margin (OPI %.3f), tier=%s", o$opi, tier)
  )
  class(out) <- "ocular_audit"
  out
}

print.ocular_audit <- function(x, ...) {
  cat("Ocular-Confound Audit\n")
  cat(sprintf("  OPI          : %+.3f   [envelope %.3f, %.3f]\n",
              x$opi, x$opi_envelope[1], x$opi_envelope[2]))
  cat(sprintf("  ceiling AUC  : %.3f   projected AUC: %.3f   p=%.3f\n",
              x$ceiling_auc, x$projected_auc, x$p_value))
  cat(sprintf("  reliability  : R2=%.3f  ->  tier = %s\n", x$reliability_r2, x$trust_tier))
  cat(sprintf("  verdict      : %s\n", x$verdict))
  invisible(x)
}

# =============================================================================
# AUDITOR RELIABILITY (Method 3)
# How well can the frozen readout reconstruct held-out gaze/ocular state in the
# TARGET montage? Reported as the mean cross-validated R^2 across ocular outputs
# (>=0 useful; <=0 means the montage carries no recoverable ocular signal).
# We retrain the readout on each target's exact channel set (montage-matched,
# NO interpolation -- Method 1) using the target's OWN ground-truth ocular
# signal where available (here: simulator / EEGEyeNet). The reliability number
# uses subject-wise CV so it never sees its own test gaze.
# =============================================================================
gaze_r2 <- function(Y_true, Y_pred) {
  Y_true <- as.matrix(Y_true); Y_pred <- as.matrix(Y_pred)
  ss_res <- colSums((Y_true - Y_pred)^2)
  ss_tot <- colSums(sweep(Y_true, 2L, colMeans(Y_true), "-")^2)
  r2 <- 1 - ss_res / ss_tot
  r2[!is.finite(r2)] <- NA_real_
  r2
}

auditor_reliability <- function(X, Y, subj, lambda = 1.0) {
  subs <- unique(subj)
  Y <- as.matrix(Y)
  pred <- matrix(NA_real_, nrow(X), ncol(Y))
  for (s in subs) {
    tr <- subj != s; te <- subj == s
    fit <- ridge_fit(X[tr, , drop = FALSE], Y[tr, , drop = FALSE], lambda)
    pred[te, ] <- ridge_predict(fit, X[te, , drop = FALSE])
  }
  r2 <- gaze_r2(Y, pred)
  list(per_output_r2 = r2, mean_r2 = mean(r2, na.rm = TRUE))
}

# Three-tier gate. Thresholds are derived empirically in run_simulation.R
# (bias-vs-fidelity sweep); the defaults here are placeholders overwritten by
# the calibrated values before any audit is interpreted.
reliability_tier <- function(mean_r2, t_conf = 0.50, t_expl = 0.20) {
  if (is.na(mean_r2))      return("inconclusive")
  if (mean_r2 >= t_conf)   return("confirmatory")
  if (mean_r2 >= t_expl)   return("exploratory")
  "inconclusive"
}

# =============================================================================
# Subject-wise (macro) AUC of a ridge cognitive decoder.
# y in {0,1}. If P is supplied, the decoder is trained & tested on the
# ocular-only projection X %*% P. AUC is computed within each held-out subject
# and averaged (macro), matching the within-subject permutation null.
# =============================================================================
loso_macro_auc <- function(X, y, subj, P = NULL, lambda = 1.0) {
  X <- as.matrix(X)
  if (!is.null(P)) X <- X %*% P
  subs <- unique(subj)
  aucs <- numeric(length(subs)); names(aucs) <- as.character(subs)
  for (i in seq_along(subs)) {
    s <- subs[i]
    tr <- subj != s; te <- subj == s
    yt <- y[te]
    if (length(unique(yt)) < 2L) { aucs[i] <- NA_real_; next }  # need both classes
    fit <- ridge_fit(X[tr, , drop = FALSE], matrix(2 * y[tr] - 1, ncol = 1), lambda)
    sc  <- as.numeric(ridge_predict(fit, X[te, , drop = FALSE]))
    aucs[i] <- auc_score(sc, yt)
  }
  mean(aucs, na.rm = TRUE)
}

# =============================================================================
# THE OCULAR PREDICTABILITY INDEX (Method 4)
# OPI = AUC(cognitive decoder on ocular-only projection) - within-subject
#       label-permutation null mean. Reported with the full-EEG AUC ceiling,
#       the null distribution, and a permutation p-value. The two-sided bias
#       margin from Method 5 is added by the caller.
# =============================================================================
compute_opi <- function(X, y, subj, P, lambda = 1.0, n_perm = 200, seed = 1) {
  set.seed(seed)
  X <- as.matrix(X)

  auc_ceiling <- loso_macro_auc(X, y, subj, P = NULL, lambda = lambda)  # full EEG
  auc_proj    <- loso_macro_auc(X, y, subj, P = P,    lambda = lambda)  # ocular only

  null <- numeric(n_perm)
  subs <- unique(subj)
  for (b in seq_len(n_perm)) {
    yp <- y
    for (s in subs) {                          # shuffle labels WITHIN subject
      idx <- which(subj == s)
      yp[idx] <- sample(y[idx])
    }
    null[b] <- loso_macro_auc(X, yp, subj, P = P, lambda = lambda)
  }
  null_mean <- mean(null, na.rm = TRUE)
  p_val <- (1 + sum(null >= auc_proj, na.rm = TRUE)) / (1 + sum(is.finite(null)))

  list(
    auc_ceiling = auc_ceiling,                 # full-EEG AUC (ceiling)
    auc_proj    = auc_proj,                    # AUC on ocular-only projection
    null_mean   = null_mean,
    null_sd     = sd(null, na.rm = TRUE),
    null_q      = quantile(null, c(.025, .5, .975), na.rm = TRUE),
    opi         = auc_proj - null_mean,        # the index
    opi_frac_of_ceiling = (auc_proj - null_mean) /
                          max(auc_ceiling - null_mean, 1e-6),
    p_value     = p_val
  )
}

# Convenience: compute the OPI using GROUND-TRUTH gaze as the "projection"
# instead of the EEG-derived ocular subspace. Used for the Method-5 bias
# envelope: the signed gap between this and compute_opi() is the measured bias.
compute_opi_from_truth <- function(gaze_features, y, subj, lambda = 1.0,
                                    n_perm = 200, seed = 1) {
  set.seed(seed)
  G <- as.matrix(gaze_features)
  auc_truth <- loso_macro_auc(G, y, subj, P = NULL, lambda = lambda)
  null <- numeric(n_perm); subs <- unique(subj)
  for (b in seq_len(n_perm)) {
    yp <- y
    for (s in subs) { idx <- which(subj == s); yp[idx] <- sample(y[idx]) }
    null[b] <- loso_macro_auc(G, yp, subj, P = NULL, lambda = lambda)
  }
  list(auc_truth = auc_truth, null_mean = mean(null, na.rm = TRUE),
       opi_truth = auc_truth - mean(null, na.rm = TRUE))
}
