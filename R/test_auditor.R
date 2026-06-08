# Unit checks for auditor.R -- run: Rscript R/test_auditor.R
suppressMessages(library(pROC))
source("R/auditor.R")
set.seed(42)
ok <- TRUE
chk <- function(name, cond) {
  cat(sprintf("[%s] %s\n", if (cond) "PASS" else "FAIL", name))
  if (!cond) ok <<- FALSE
}

# 1. AUC matches pROC (within numerical tolerance), incl. ties
s <- rnorm(500); l <- rbinom(500, 1, plogis(0.8 * s))
a_mine <- auc_score(s, l)
a_proc <- as.numeric(pROC::auc(pROC::roc(l, s, quiet = TRUE, direction = "<")))
chk(sprintf("AUC vs pROC (%.4f vs %.4f)", a_mine, a_proc), abs(a_mine - a_proc) < 1e-6)
st <- sample(1:5, 500, TRUE)                       # heavy ties
at_mine <- auc_score(st, l)
at_proc <- as.numeric(pROC::auc(pROC::roc(l, st, quiet = TRUE, direction = "<")))
chk(sprintf("AUC with ties vs pROC (%.4f vs %.4f)", at_mine, at_proc), abs(at_mine - at_proc) < 1e-6)

# 2. Projector is symmetric & idempotent, with the right rank
W <- matrix(rnorm(20 * 3), 20, 3)                  # 20 features, 3 ocular outputs
pj <- ocular_projection(W)
P <- pj$P
chk("Projector symmetric", max(abs(P - t(P))) < 1e-10)
chk("Projector idempotent (P^2 = P)", max(abs(P %*% P - P)) < 1e-9)
chk(sprintf("Projector rank = 3 (got %d)", pj$rank), pj$rank == 3)
chk("Projector fixes its own column space", {
  v <- W %*% rnorm(3); max(abs(P %*% v - v)) < 1e-9
})
chk("Projector kills the orthogonal complement", {
  # a vector orthogonal to col(W) must map to ~0
  qrW <- qr.Q(qr(W)); comp <- diag(20) - qrW %*% t(qrW)
  v <- comp %*% rnorm(20); max(abs(P %*% v)) < 1e-9
})

# 3. Ridge recovers a known linear map at small lambda
Xs <- matrix(rnorm(800 * 6), 800, 6); Btrue <- matrix(rnorm(6 * 2), 6, 2)
Ys <- Xs %*% Btrue + matrix(rnorm(800 * 2, sd = 0.01), 800, 2)
fit <- ridge_fit(Xs, Ys, lambda = 1e-4)
chk(sprintf("Ridge recovers weights (max err %.3f)", max(abs(fit$W - Btrue))),
    max(abs(fit$W - Btrue)) < 0.05)

# 4. End-to-end: when the cognitive signal lives INSIDE the ocular subspace,
#    OPI should be clearly positive; when it lives OUTSIDE, OPI ~ 0.
n_sub <- 6; per <- 80; N <- n_sub * per
subj <- rep(1:n_sub, each = per); y <- rbinom(N, 1, 0.5)
# ocular outputs read channels 1:3; cognitive signal A is in ch1, signal B in ch10
Xbase <- matrix(rnorm(N * 12), N, 12)
Wsrc <- matrix(0, 12, 3); Wsrc[1,1] <- 1; Wsrc[2,2] <- 1; Wsrc[3,3] <- 1
P4 <- ocular_projection(Wsrc)$P
Xin  <- Xbase; Xin[,1]  <- Xin[,1]  + 1.2 * (y - 0.5)   # signal inside subspace
Xout <- Xbase; Xout[,10] <- Xout[,10] + 1.2 * (y - 0.5)  # signal outside subspace
opi_in  <- compute_opi(Xin,  y, subj, P4, n_perm = 100, seed = 1)
opi_out <- compute_opi(Xout, y, subj, P4, n_perm = 100, seed = 1)
cat(sprintf("    OPI(signal inside subspace)  = %+.3f  (proj AUC %.3f, p=%.3f)\n",
            opi_in$opi, opi_in$auc_proj, opi_in$p_value))
cat(sprintf("    OPI(signal outside subspace) = %+.3f  (proj AUC %.3f, p=%.3f)\n",
            opi_out$opi, opi_out$auc_proj, opi_out$p_value))
chk("OPI positive when signal is in the ocular subspace", opi_in$opi > 0.10)
# Outside-subspace OPI must be non-significant and far below the inside case.
# (Sign can be slightly negative from finite-sample macro-AUC noise; what
#  matters is p ~ NS and a large gap from the inside-subspace OPI.)
chk("OPI not significant when signal is outside the ocular subspace",
    opi_out$p_value > 0.20)
chk("OPI inside >> OPI outside", (opi_in$opi - opi_out$opi) > 0.20)

cat(if (ok) "\nALL TESTS PASSED\n" else "\nSOME TESTS FAILED\n")
quit(status = if (ok) 0 else 1)
