# demo_tool.R -- minimal usage of the shipped auditor (Method 7).
# Train a frozen auditor on a gaze-rich source, then audit two target datasets.
source("R/auditor.R"); source("R/simulate.R")

# 1. Train & freeze the auditor (subject= -> stores montage reliability).
src <- simulate_dataset("neutral", n_subjects = 20, trials_per = 200,
                        sacc_deg = 8, montage = "dense", seed = 100)
auditor <- train_auditor(src$X, src$gaze, lambda = 1, subject = src$subj)
cat(sprintf("Frozen auditor: rank %d, montage reliability R2 = %.3f\n\n",
            auditor$rank, auditor$reliability_r2))

# 2. Audit a CONTAMINATED dataset (gaze coupled to the cognitive label).
bad <- simulate_dataset("ocular_corr", coupling = 0.8, sacc_deg = 8,
                        montage = "dense", seed = 7)
cat(">> Dataset A (eyes coupled to condition):\n")
print(audit_dataset(bad$X, bad$y, bad$subj, auditor, n_perm = 200))

# 3. Audit a CLEAN dataset (covert attention; neural signal, no ocular coupling).
good <- simulate_dataset("covert", sacc_deg = 8, montage = "dense", seed = 7)
cat("\n>> Dataset B (covert attention, no ocular coupling):\n")
print(audit_dataset(good$X, good$y, good$subj, auditor, n_perm = 200))
