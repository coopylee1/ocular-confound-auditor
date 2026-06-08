# =============================================================================
# run_leakage.R  --  E6: the neural-leakage term is governed by the principal
# angle between the neural cognitive topography and the ocular readout span.
#
# Tests the proposed bias-decomposition claim:
#   OPI = (true ocular-decodable component)  -  (under-capture deficit)
#                                            +  (neural-leakage term)
# where the neural-leakage term is BOUNDED by the subspace overlap
#   s = ||P t|| / ||t||  =  cos(principal angle) between the neural contrast
# topography t and the frozen ocular subspace P. If true, leakage -> 0 as the
# neural signal becomes orthogonal to the ocular subspace (posterior sources),
# and grows as it aligns (frontal sources). No ocular coupling is present, so any
# positive OPI here is PURE neural leakage.
# =============================================================================
suppressMessages(library(parallel))
source("R/auditor.R"); source("R/simulate.R")
set.seed(1)

mont <- make_montage(); C <- mont$coords[mont$full, , drop = FALSE]
topo_all <- make_topographies(list(coords = C, eye_left = mont$eye_left,
                                   eye_right = mont$eye_right))
# Frozen ocular auditor (trained on gaze-rich neutral source, as elsewhere)
src <- simulate_dataset("neutral", n_subjects = 20, trials_per = 200,
                        sacc_deg = 8, montage = "full", seed = 100)
AUD <- train_auditor(src$X, src$gaze, lambda = 1)
P <- AUD$P

# A lateralised neural contrast source at vertical position ypos (front=+, back=-)
neural_contrast_topo <- function(ypos) {
  tr <- .prox(C, c( 0.50, ypos), 0.28)
  tl <- .prox(C, c(-0.50, ypos), 0.28)
  tr - tl                                   # right-minus-left contrast direction
}
overlap <- function(t) sqrt(sum((P %*% t)^2)) / sqrt(sum(t^2))   # cos(principal angle)

gen_neural_only <- function(ypos, seed, n_subjects = 10, trials_per = 120,
                            neural_uV = 2.2) {
  set.seed(seed)
  N <- n_subjects * trials_per; subj <- rep(seq_len(n_subjects), each = trials_per)
  y <- rbinom(N, 1, 0.5); side <- ifelse(y == 1, 1, -1)
  t <- neural_contrast_topo(ypos)
  X <- .bg_noise(N, nrow(C), C, sd_smooth = 3, sd_white = 3)
  X <- X + outer(neural_uV * side, t)                  # neural signal, label-lateralised
  # independent ocular activity NOT coupled to the label (so OPI=pure leakage)
  gx <- rnorm(N, 0, 8) * sample(c(-1,1), N, TRUE)
  X <- X + outer(gx * 16, topo_all$heog) + outer(rnorm(N,0,2.5)*6, topo_all$veog)
  list(X = X, y = y, subj = subj)
}

ypositions <- seq(0.55, -0.90, length.out = 10)   # frontal (aligned) -> posterior (orthogonal)
seeds <- 1:8
rows <- mclapply(ypositions, function(yp) {
  ov <- overlap(neural_contrast_topo(yp))
  opis <- sapply(seeds, function(s) {
    d <- gen_neural_only(yp, seed = s * 101 + 5)
    compute_opi(d$X, d$y, d$subj, P, n_perm = 150, seed = s)$opi
  })
  data.frame(ypos = yp, overlap = ov, opi_leak = mean(opis), opi_sd = sd(opis))
}, mc.cores = max(1, min(detectCores()-1, 8)))
E6 <- do.call(rbind, rows)
E6 <- E6[order(E6$overlap), ]
cat("=== E6: neural-leakage OPI vs ocular-subspace overlap (principal angle) ===\n")
print(E6, digits = 3, row.names = FALSE)
cat(sprintf("\nSpearman rho(overlap, leakage) = %.3f  (expect ~+1: leakage tracks overlap)\n",
            cor(E6$overlap, E6$opi_leak, method = "spearman")))
cat(sprintf("leakage at min overlap (%.2f) = %+.3f ;  at max overlap (%.2f) = %+.3f\n",
            min(E6$overlap), E6$opi_leak[which.min(E6$overlap)],
            max(E6$overlap), E6$opi_leak[which.max(E6$overlap)]))
saveRDS(E6, "results/tables/E6_leakage.rds")
write.csv(E6, "results/tables/E6_leakage.csv", row.names = FALSE)
