# =============================================================================
# simulate.R  --  Grounded EEG + gaze forward-model simulator
#
# Why a simulator (and why it is NOT a cop-out): the paper's own calibration
# instruments (Methods 5-6, H3) are built on SYNTHETIC condition-correlated
# ocular injection -- it is the only way to obtain ground-truth knowledge of how
# much of a decodable signal is ocular. The simulator gives us that ground
# truth so we can prove the auditor recovers it. Every amplitude/topography
# below is grounded in the fact-checked literature:
#
#   * Corneo-retinal standing potential ~0.4-1.0 mV at the eye; scalp HEOG
#     deflection ~ 16 uV per degree of horizontal rotation near the canthi
#     (Malmivuo & Plonsey, Bioelectromagnetism ch.28; standard EOG figure).
#   * Saccadic spike potential: brief (~22 ms) biphasic transient at saccade
#     onset, broadband ~20-90 Hz, focal at periocular sites, lateralised to the
#     saccade target, amplitude GROWS with saccade size (Keren, Yuval-Greenberg
#     & Deouell 2010, NeuroImage 49:2248). We model it as a few-uV periocular
#     source scaling with saccade amplitude.
#   * Blink: large (~50-100 uV) frontopolar monophasic deflection.
#   * N2pc: lateralised posterior negativity ~1-3 uV, contralateral to the
#     attended hemifield, peaking at PO7/PO8 (Luck; ERP CORE).
#
# The demonstrator works on SPATIAL (channel-space) trial features -- one
# amplitude per channel per trial. That is exactly the object the LOCKED
# operator (an orthogonal projection in channel space) acts on. The full
# instrument would use spatiotemporal windows; the spatial subspace is its
# essence and is what we validate here. This simplification is stated in README.
# =============================================================================

# ---- Montage: full 64-channel 10-10 layout (x = right+, y = front+) -----------
# The auditor operates on ALL available channels. We model a standard 64-channel
# 10-10 cap (the "full" montage). The 30 ERP-CORE channels are an explicitly
# labelled SUBSET of these 64 (used only to characterise how channel count /
# near-eye frontal coverage affects the auditor -- i.e. the reliability gate). We
# never silently reduce: subsetting is a deliberately studied variable, and the
# real-data audit (eegmmidb) likewise uses all 64 channels.
make_montage <- function() {
  erpcore <- rbind(                                   # the 30-channel ERP-CORE subset
    Fp1=c(-0.30,0.95), Fp2=c(0.30,0.95),
    F7=c(-0.75,0.60), F3=c(-0.35,0.62), Fz=c(0,0.65), F4=c(0.35,0.62), F8=c(0.75,0.60),
    FC3=c(-0.40,0.35), FCz=c(0,0.38), FC4=c(0.40,0.35),
    C5=c(-0.70,0.05), C3=c(-0.38,0.05), Cz=c(0,0.05), C4=c(0.38,0.05), C6=c(0.70,0.05),
    CPz=c(0,-0.25),
    P7=c(-0.75,-0.50), P3=c(-0.38,-0.50), Pz=c(0,-0.50), P4=c(0.38,-0.50), P8=c(0.75,-0.50),
    P9=c(-0.80,-0.62), P10=c(0.80,-0.62),
    PO7=c(-0.55,-0.72), PO3=c(-0.30,-0.75), PO4=c(0.30,-0.75), PO8=c(0.55,-0.72),
    O1=c(-0.25,-0.95), Oz=c(0,-0.98), O2=c(0.25,-0.95)
  )
  extra <- rbind(                                     # remaining 34 channels of the 64-ch cap
    Fpz=c(0,0.97),
    AF7=c(-0.55,0.80), AF3=c(-0.28,0.83), AFz=c(0,0.84), AF4=c(0.28,0.83), AF8=c(0.55,0.80),
    F5=c(-0.55,0.62), F1=c(-0.18,0.67), F2=c(0.18,0.67), F6=c(0.55,0.62),
    FT7=c(-0.78,0.33), FC5=c(-0.58,0.36), FC1=c(-0.19,0.40), FC2=c(0.19,0.40),
    FC6=c(0.58,0.36), FT8=c(0.78,0.33),
    T7=c(-0.82,0.05), C1=c(-0.20,0.05), C2=c(0.20,0.05), T8=c(0.82,0.05),
    TP7=c(-0.78,-0.27), CP5=c(-0.58,-0.25), CP3=c(-0.38,-0.24), CP1=c(-0.19,-0.24),
    CP2=c(0.19,-0.24), CP4=c(0.38,-0.24), CP6=c(0.58,-0.25), TP8=c(0.78,-0.27),
    P5=c(-0.55,-0.52), P1=c(-0.18,-0.52), P2=c(0.18,-0.52), P6=c(0.55,-0.52),
    POz=c(0,-0.77), Iz=c(0,-1.05)
  )
  all <- rbind(erpcore, extra)
  colnames(all) <- c("x","y")
  list(coords = all,
       erpcore = rownames(erpcore),                   # 30-channel subset
       dense   = rownames(all),                        # alias: the full 64-channel montage
       full    = rownames(all),
       eye_left = c(-0.30,1.15), eye_right = c(0.30,1.15))
}

.prox <- function(coords, point, sigma = 0.45) {
  exp(-((coords[,"x"]-point[1])^2 + (coords[,"y"]-point[2])^2) / (2*sigma^2))
}

# Spatial topographies (length = n channels) for each ocular / neural source.
make_topographies <- function(mont) {
  C  <- mont$coords; el <- mont$eye_left; er <- mont$eye_right
  pl <- .prox(C, el, 0.40); pr <- .prox(C, er, 0.40)            # near each eye
  list(
    heog  = (pr - pl),                                          # horizontal: antisymmetric
    veog  = (pl + pr),                                          # vertical / blink: frontopolar symmetric
    sp_l  = .prox(C, el, 0.22),                                 # saccadic spike, focal at left eye
    sp_r  = .prox(C, er, 0.22),                                 # focal at right eye
    blink = .prox(C, c(0,1.15), 0.30),                          # blink: focal mid-frontopolar
    n2pc_l= .prox(C, c(-0.55,-0.72), 0.28),                     # left posterior (PO7)
    n2pc_r= .prox(C, c( 0.55,-0.72), 0.28),                     # right posterior (PO8)
    muC3  = .prox(C, c(-0.38,0.05), 0.20),                      # left sensorimotor (focal mu/beta ERD)
    muC4  = .prox(C, c( 0.38,0.05), 0.20)                       # right sensorimotor (focal mu/beta ERD)
  )
}

# Spatially smooth background noise: shared low-rank field + per-channel white.
.bg_noise <- function(n_trials, K, coords, sd_smooth = 3.0, sd_white = 3.0) {
  # 4 smooth spatial modes (random gaussian bumps) shared structure across chans
  centers <- cbind(runif(4,-0.8,0.8), runif(4,-0.9,0.9))
  B <- sapply(1:4, function(i) .prox(coords, centers[i,], 0.5))  # K x 4
  amp <- matrix(rnorm(n_trials*4, sd = sd_smooth), n_trials, 4)
  amp %*% t(B) + matrix(rnorm(n_trials*K, sd = sd_white), n_trials, K)
}

# =============================================================================
# simulate_dataset(): one scenario.
#   scenario:
#     "ocular_corr"  - gaze direction coupled to label (condition-correlated
#                      ocular)  ==> positive control, OPI should be HIGH
#     "covert"       - neural N2pc lateralised by label, gaze CENTRAL (no ocular
#                      coupling)  ==> tests neural leakage, OPI should be ~NULL
#     "motor"        - non-spatial label, central mu-rhythm neural source, gaze
#                      independent  ==> negative control, OPI ~NULL
#     "blink_corr"   - DIFFERENT ARTIFACT: blink rate coupled to the label, gaze
#                      independent  ==> positive control for blink contamination
#     "neutral"      - neural posterior signal + independent random eye movements
#                      ==> baseline negative control
#   coupling   in [0,1]: strength of gaze<->label (ocular_corr) or blink<->label coupling
#   sacc_deg   typical saccade amplitude in degrees (amplitude sweep / regime)
#   montage    "full" (64 channels, default) or "erpcore" (30-channel subset)
# Returns EEG trial-features X (trials x channels), labels y, subject ids, and
# the ground-truth ocular targets (for training the readout / bias envelope).
# =============================================================================
simulate_dataset <- function(scenario = "ocular_corr",
                              n_subjects = 12, trials_per = 120,
                              coupling = 0.8, sacc_deg = 8, blink_rate = 0.12,
                              montage = "dense",
                              neural_uV = 1.8, heog_uV_per_deg = 16,
                              sp_uV_per_deg = 0.7, blink_uV = 80,
                              noise_white = 3.0, noise_smooth = 3.0,
                              seed = 1) {
  set.seed(seed)
  mont <- make_montage()
  use  <- if (montage == "erpcore") mont$erpcore else mont$full   # "full"/"dense" -> 64 ch
  C    <- mont$coords[use, , drop = FALSE]
  topo <- make_topographies(list(coords = C, eye_left = mont$eye_left,
                                 eye_right = mont$eye_right))
  K <- length(use)
  N <- n_subjects * trials_per
  subj <- rep(seq_len(n_subjects), each = trials_per)

  y    <- rbinom(N, 1, 0.5)                     # cognitive label (balanced)
  side <- ifelse(y == 1, 1, -1)                 # attended/target hemifield

  # ---- Gaze behaviour --------------------------------------------------------
  gx <- numeric(N); gy <- numeric(N)
  sacc_amp <- numeric(N); sacc_dir <- numeric(N)
  if (scenario == "ocular_corr") {
    # gaze deviates toward the labelled side with prob increasing in coupling
    p_toward <- 0.5 + 0.5 * coupling
    toward   <- rbinom(N, 1, p_toward) * 2 - 1            # +1 follow side, -1 against
    dir      <- side * toward
    amp      <- abs(rnorm(N, sacc_deg, sacc_deg * 0.25))
    gx       <- dir * amp
    # saccade amplitude correlates with, but is not identical to, gaze position
    sacc_amp <- amp + abs(rnorm(N, 0, sacc_deg * 0.20)); sacc_dir <- dir
  } else if (scenario == "covert") {
    # covert attention: eyes stay near centre regardless of label (small noise)
    gx <- rnorm(N, 0, 0.3); sacc_amp <- abs(rnorm(N, 0.3, 0.15))
    sacc_dir <- sign(rnorm(N))
  } else if (scenario == "motor") {
    gx <- rnorm(N, 0, sacc_deg * 0.4)
    sacc_amp <- abs(gx) + abs(rnorm(N, 0, sacc_deg * 0.2)); sacc_dir <- sign(gx)
  } else if (scenario == "blink_corr") {
    # DIFFERENT ARTIFACT TYPE: eyes roughly central & independent of the label,
    # but BLINK rate is coupled to the label (handled in the blink block below).
    gx <- rnorm(N, 0, 0.4); sacc_amp <- abs(rnorm(N, 0.3, 0.15)); sacc_dir <- sign(rnorm(N))
  } else { # neutral
    gx <- rnorm(N, 0, sacc_deg) * sample(c(-1,1), N, TRUE)
    sacc_amp <- abs(gx) + abs(rnorm(N, 0, sacc_deg * 0.2)); sacc_dir <- sign(gx)
  }
  gy <- rnorm(N, 0, 2.5)                          # vertical gaze: a real, separable dimension
  blink_p <- rep(blink_rate, N)                   # per-trial blink probability
  if (scenario == "blink_corr")                   # condition-correlated blink artifact
    blink_p <- pmin(0.90, blink_rate + coupling * 0.45 * y)
  blink <- rbinom(N, 1, blink_p)

  # ---- Assemble EEG = sum of source topographies * per-trial amplitudes ------
  X <- .bg_noise(N, K, C, sd_smooth = noise_smooth, sd_white = noise_white)

  # ocular: HEOG (corneo-retinal), VEOG, saccadic spike, blink
  X <- X + outer(gx * heog_uV_per_deg, topo$heog)
  X <- X + outer(gy * 6,               topo$veog)
  sp_topo <- ifelse(sacc_dir >= 0, 1, 0)                  # toward-right uses right-eye spike
  X <- X + outer(sacc_amp * sp_uV_per_deg * (sacc_dir>0), topo$sp_r) +
           outer(sacc_amp * sp_uV_per_deg * (sacc_dir<0), topo$sp_l)
  X <- X + outer(blink * blink_uV, topo$blink)

  # neural cognitive source
  if (scenario == "motor") {
    # label encoded by lateralised mu desync: y=1 -> stronger right-central
    X <- X + outer(neural_uV * (y - 0.5) * 2, topo$muC4) -
             outer(neural_uV * (y - 0.5) * 2, topo$muC3)
  } else {
    # N2pc: negativity contralateral to attended side
    contra <- ifelse(side > 0, 1, 0)                      # side=right -> left PO neg
    X <- X - outer(neural_uV * (side > 0), topo$n2pc_l) -
             outer(neural_uV * (side < 0), topo$n2pc_r)
  }

  colnames(X) <- use
  # Ground-truth ocular targets the auditor's readout is trained to reconstruct
  gaze_truth <- cbind(gx = gx, gy = gy,
                      sacc = sacc_amp * sacc_dir,
                      blink = blink)
  list(X = X, y = y, subj = subj, side = side, gaze = gaze_truth,
       channels = use, scenario = scenario, coupling = coupling,
       sacc_deg = sacc_deg, montage = montage)
}
