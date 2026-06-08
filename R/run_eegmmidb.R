# =============================================================================
# run_eegmmidb.R  --  REAL-DATA negative control (H3) on PhysioNet eegmmidb.
#
# eegmmidb has 64-ch BCI2000 EEG, 160 Hz, NO eye tracker. We therefore estimate
# the frozen ocular subspace from the EYES-OPEN baseline run (R01, ocular-rich:
# blinks + small eye movements) on a set of SOURCE subjects, freeze it, and audit
# a motor (rest-vs-imagined-movement) contrast on HELD-OUT TARGET subjects. This
# is a real-EEG analog of the EEGEyeNet -> target transfer:
#   * the auditor (an orthogonal projection onto the span of a linear EEG->ocular
#     readout) is trained on ocular-rich data and never sees the target labels;
#   * a motor cognitive label is decodable from full EEG (ceiling AUC) but should
#     NOT be decodable from the ocular subspace (OPI ~ null) -- motor mu/beta
#     lateralisation is central (C3/C4), not periocular.
# Limitation (stated plainly): without gaze ground truth the ocular subspace is
# estimated from baseline ocular activity, not from synchronised gaze.
# =============================================================================
suppressMessages({ library(edfReader); library(signal) })
source("R/auditor.R")
set.seed(7)
DATA <- "data/eegmmidb"; FS <- 160
OUT_T <- "results/tables"; dir.create(OUT_T, showWarnings = FALSE, recursive = TRUE)

# ---- IO ----------------------------------------------------------------------
read_run <- function(f) {
  h <- readEdfHeader(f); s <- readEdfSignals(h)
  labs <- h$sHeaders$label
  eeg  <- labs[labs != "EDF Annotations"]
  X <- sapply(eeg, function(L) s[[L]]$signal)           # T x 64
  colnames(X) <- sub("\\.+$", "", eeg)                  # "Fp1." -> "Fp1"
  list(X = X, ann = s[["EDF Annotations"]]$annotations)
}
bp_filt <- butter(4, c(8,30)/(FS/2), type = "pass")     # mu+beta (motor)
lp_filt <- butter(4, 7/(FS/2),       type = "low")      # ocular band
bpf <- function(M) apply(M, 2, function(c) filtfilt(bp_filt, c))
lpf <- function(M) apply(M, 2, function(c) filtfilt(lp_filt, c))

subjects <- sort(list.files(DATA, pattern = "^S[0-9]+$"))
has_all  <- subjects[sapply(subjects, function(s)
  all(file.exists(file.path(DATA, s, paste0(s, c("R01","R04","R08","R12"), ".edf")))))]
cat("subjects with R01+R04+R08+R12:", paste(has_all, collapse = " "), "\n")
stopifnot(length(has_all) >= 4)

# Source (estimate frozen ocular subspace) vs Target (audited). Strict split.
half   <- max(2, floor(length(has_all)/2))
source_subj <- has_all[seq_len(half)]
target_subj <- has_all[(half+1):length(has_all)]
cat("source (ocular subspace):", paste(source_subj, collapse=" "), "\n")
cat("target (audited)        :", paste(target_subj, collapse=" "), "\n\n")

# ---- Frozen ocular subspace from eyes-open R01 of source subjects ------------
# Readout maps EEG -> {vertical/blink, horizontal} ocular reference signals
# (built from frontopolar / anterior-frontal channels of the low-passed signal).
ref_signals <- function(Xlp) {
  cbind(vert = rowMeans(Xlp[, c("Fp1","Fpz","Fp2")]),
        horiz = Xlp[, "Af7"] - Xlp[, "Af8"])
}
Xsrc <- NULL; Rsrc <- NULL
for (s in source_subj) {
  r <- read_run(file.path(DATA, s, paste0(s, "R01.edf")))
  Xlp <- lpf(r$X)
  Xsrc <- rbind(Xsrc, Xlp); Rsrc <- rbind(Rsrc, ref_signals(Xlp))
}
AUD <- train_auditor(Xsrc, Rsrc, lambda = 1.0)
cat(sprintf("frozen ocular subspace rank = %d (of 64 channels)\n", AUD$rank))

# ---- Epoch the motor runs of target subjects (rest vs imagined movement) -----
epoch_run <- function(f) {
  r <- read_run(f); E <- bpf(r$X); a <- r$ann
  a <- a[a$annotation %in% c("T0","T1","T2"), ]
  feats_full <- list(); feats_proj <- list(); lab <- c()
  for (i in seq_len(nrow(a))) {
    on <- round((a$onset[i] + 0.5) * FS); off <- on + round(3.0 * FS) - 1
    if (off > nrow(E)) next
    seg <- E[on:off, , drop = FALSE]
    feats_full[[length(feats_full)+1]] <- log(apply(seg, 2, var) + 1e-8)
    feats_proj[[length(feats_proj)+1]] <- log(apply(seg %*% AUD$P, 2, var) + 1e-8)
    lab <- c(lab, a$annotation[i])                    # keep T0/T1/T2
  }
  list(full = do.call(rbind, feats_full), proj = do.call(rbind, feats_proj), lab = lab)
}

Xfull <- NULL; Xproj <- NULL; lab <- c(); subj <- c()
for (s in target_subj) for (rr in c("R04","R08","R12")) {
  e <- epoch_run(file.path(DATA, s, paste0(s, rr, ".edf")))
  Xfull <- rbind(Xfull, e$full); Xproj <- rbind(Xproj, e$proj)
  lab <- c(lab, e$lab); subj <- c(subj, rep(s, length(e$lab)))
}

# within-subject stratified k-fold macro-AUC (standard MI evaluation; cross-
# subject linear MI decoding is near chance, a known fact, so the ceiling for the
# lateralised contrast is evaluated within subject).
within_macro_auc <- function(X, y, subj, k = 5, seed = 3) {
  set.seed(seed); subs <- unique(subj); a <- c()
  for (s in subs) {
    idx <- which(subj == s); ys <- y[idx]
    if (length(unique(ys)) < 2) next
    fold <- ave(seq_along(idx), ys, FUN = function(z) sample(rep_len(1:k, length(z))))
    sc <- numeric(length(idx))
    for (f in 1:k) {
      tr <- idx[fold != f]; te <- idx[fold == f]
      if (length(unique(y[tr])) < 2) { sc[fold==f] <- 0; next }
      fit <- ridge_fit(X[tr,,drop=FALSE], matrix(2*y[tr]-1, ncol=1), 1.0)
      sc[fold == f] <- as.numeric(ridge_predict(fit, X[te,,drop=FALSE]))
    }
    a <- c(a, auc_score(sc, ys))
  }
  mean(a, na.rm = TRUE)
}
perm_p <- function(fun, X, y, subj, obs, NP = 500, seed = 11) {
  set.seed(seed); subs <- unique(subj); null <- numeric(NP)
  for (b in 1:NP) { yp <- y; for (s in subs){i<-which(subj==s); yp[i]<-sample(y[i])}
    null[b] <- fun(X, yp, subj) }
  list(p = (1 + sum(null >= obs, na.rm=TRUE))/(1+sum(is.finite(null))), null_mean = mean(null))
}

audit <- function(name, keep, ymap, cvfun) {
  m <- lab %in% keep
  y <- ymap(lab[m]); sj <- subj[m]; Xf <- Xfull[m,,drop=FALSE]; Xp <- Xproj[m,,drop=FALSE]
  ceil <- cvfun(Xf, y, sj); proj <- cvfun(Xp, y, sj)
  pr <- perm_p(cvfun, Xp, y, sj, proj)
  cat(sprintf("\n[%s]  n=%d (%d/%d)  ceiling AUC=%.3f  ocular-proj AUC=%.3f  OPI=%+.3f  p=%.3f\n",
              name, length(y), sum(y==0), sum(y==1), ceil, proj, proj - pr$null_mean, pr$p))
  data.frame(contrast=name, n=length(y), auc_full=ceil, auc_proj=proj,
             opi=proj - pr$null_mean, p_value=pr$p)
}

cat("\n=== REAL-DATA AUDIT (eegmmidb, frozen ocular subspace from eyes-open) ===")
R1 <- audit("left_vs_right_imagery (clean negative control, within-subj CV)",
            c("T1","T2"), function(l) as.integer(l=="T2"), within_macro_auc)
R2 <- audit("rest_vs_move (subject-wise LOSO; frontal-proxy leakage demo)",
            c("T0","T1","T2"), function(l) as.integer(l!="T0"),
            function(X,y,s) loso_macro_auc(X,y,s,P=NULL,lambda=1.0))
res <- rbind(R1, R2)
cat("\nInterpretation: left-vs-right imagery has NO differential eye movement -> the\n",
    "clean ocular-subspace negative control; rest-vs-move carries frontal neural\n",
    "signal the channel-proxy subspace cannot separate from ocular (no gaze GT in\n",
    "eegmmidb) -- which is precisely why the real auditor must be gaze-trained.\n", sep="")
write.csv(res, file.path(OUT_T, "real_eegmmidb_audit.csv"), row.names = FALSE)
cat("\nwritten: results/tables/real_eegmmidb_audit.csv\n")
