# Quantifying Ocular Confounds in EEG Decoding Without an Eye Tracker — execution report

> **Superseded — pre-revision execution record.** This report documents the original run and
> uses the paper's earlier framing ("ocular subspace", "calibrated", "ceiling", a three-tier
> gate). The peer-reviewed claims are narrower; see
> [`camera_ready/manuscript.pdf`](camera_ready/manuscript.pdf) and the E7–E10 tables in
> `results/tables/` for the baselines, per-output reliability, and held-out calibration added
> in revision.


**A portable ocular-confound auditor for EEG decoding: implementation + in-silico
and real-data validation.**

This report documents what was actually built and run for Paper 15. Every number
below was produced by executing the R code in `R/`; nothing is asserted without a
run behind it. Source of truth: `results/tables/*.csv`, `results/figures/*.png`.

Date run: 2026-06-05 · R 4.6.0 · 12-core CPU. Reproduce with `bash run_all.sh`.

---

## 1. What this instrument is (and why it is faithfully implementable in R)

The proposal's **locked operator** (Method 1) is *"the orthogonal projection onto
the span of the frozen linear gaze/ocular-state readout weights."* The **index**
(Method 4) is *"the AUC of a paradigm-matched cognitive decoder applied to that
ocular-only projection, minus a within-subject label-permutation null."* Both are
linear algebra + permutation statistics. A deep EEGEyeNet network is one way to
*train* a gaze readout, but the object that is **frozen and shipped** is the linear
readout and its projector — so an R implementation built on a ridge readout is not
a shortcut, it *is* the locked instrument. The auditor is therefore implemented in
full in [`R/auditor.R`](R/auditor.R):

- `train_auditor()` — ridge readout EEG→{gaze_x, gaze_y, saccade, blink}, then
  `ocular_projection()` = `U Uᵀ`, the orthogonal projector onto the column space of
  the readout weights (`U` = its left singular vectors). **The projection operator
  is locked to this single definition before any audit.**
- `compute_opi()` — subject-wise (macro) AUC of a ridge cognitive decoder on the
  ocular-only projection `X·P`, minus the within-subject label-permutation null,
  with the **full-EEG AUC as ceiling** and a permutation p-value.
- `auditor_reliability()` / `reliability_tier()` — cross-validated gaze-reconstruction
  R² → confirmatory / exploratory / inconclusive gate.
- `compute_opi_from_truth()` — the same index from **ground-truth gaze**, used to
  measure the bias envelope (Method 5).

The core math is unit-tested in [`R/test_auditor.R`](R/test_auditor.R): AUC matches
`pROC` to 1e-6 (incl. ties), the projector is symmetric/idempotent/correct-rank,
ridge recovers known weights, and the OPI is large when the cognitive signal lives
in the ocular subspace and non-significant when it does not. **All tests pass.**

## 2. The simulator (grounded, not arbitrary)

A synthetic EEG+gaze forward model ([`R/simulate.R`](R/simulate.R)) is the only way
to get the **ground-truth knowledge of "how much is ocular"** that calibration
needs — and synthetic condition-correlated injection is itself one of the proposal's
instruments (Methods 5–6). Every amplitude/topography is grounded in the
fact-checked literature: corneo-retinal HEOG ≈ 16 µV/deg, Keren-characterized
saccadic spike (focal periocular, scales with saccade size), ~80 µV frontopolar
blink, ~1.8 µV lateralized posterior N2pc. The montage is the exact 30-channel
ERP-CORE set, optionally extended with 8 periocular/frontopolar channels that a
128-ch EGI net has and a 30-ch cap lacks (Figure 1). Source topographies are shown
in Figure 2.

## 3. Results

### H3 — graded sensitivity (E1 battery, E2 sweep) ✅

Frozen auditor trained on a large-saccade source; **full 64-channel montage**; 12
independent simulated "studies" per cell; 200-permutation null.
`results/tables/E1_battery_summary.csv`:

| scenario | mean OPI | proj AUC | ceiling AUC | frac p<.05 |
|---|---|---|---|---|
| **saccade-coupled ocular** (positive control) | **+0.429** | 0.929 | 0.930 | **1.00** |
| **blink-coupled ocular** (positive control, distinct artifact) | **+0.198** | 0.697 | 0.781 | **1.00** |
| motor imagery | +0.031 | 0.531 | 0.842 | 0.42 |
| covert attention | +0.007 | 0.507 | 0.785 | 0.17 |
| label-shuffle | +0.007 | 0.507 | 0.505 | 0.00 |
| neutral | −0.003 | 0.496 | 0.635 | 0.00 |

Two positive controls test **distinct artifact types**: saccade-coupled (0.429) and
**blink-coupled (0.198)** ocular contamination are both detected at 100% significance —
the auditor catches blinks as well as saccades. The **covert-attention** control is
decisive: the label is decodable from full EEG (ceiling 0.79) but **not** from the
ocular projection (OPI 0.007, NS) — the auditor does **not** absorb lateralized posterior
neural signal. The four negative controls (≤0.031) set an empirical **calibrated
neural-leakage margin ≈ 0.05** (Figure 3) that the index must clear, beyond the bare
permutation null. (Modelling motor imagery as a *focal* sensorimotor source keeps its
overlap with the ocular subspace small; a broader lateralised source would leak more, as
the principal-angle analysis below quantifies.)

The coupling sweep (`E2_coupling_summary.csv`, Figure 4) shows the index is a *meter*,
rising monotonically and near-linearly with gaze↔label coupling ρ:

| ρ | 0.0 | 0.2 | 0.4 | 0.6 | 0.8 | 1.0 |
|---|---|---|---|---|---|---|
| OPI | −0.014 | 0.100 | 0.202 | 0.303 | 0.404 | 0.500 |

### Method 5 — measured two-sided bias envelope (E3) ✅

Across the saccade-amplitude range, the auditor's OPI is compared to the OPI computed
from **ground-truth gaze** (Figure 5, `E3_amplitude_summary.csv`). From 0.5° upward the
bias is within ±0.01 on both montages; only at the **fixational extreme (0.25°)** does it
grow — to **−0.034** on the full 64-channel montage and **−0.056** on the 30-channel
subset — i.e. the auditor **under-states** the ocular contribution (the *conservative*
direction). No over-statement was observed anywhere. The "surrogate ≤ true" guarantee
thus holds in the regime where it matters, and the magnitude of the breach is
**measured**, not assumed — and folded into the gate as a two-sided margin.

### Method 3 — empirically calibrated reliability gate (E4) ✅

From the bias-vs-fidelity relationship we derive the thresholds (not hand-set):
**t_conf = 0.664, t_expl = 0.497** (mean gaze R²; `calibration.json`). Tiering the
amplitude×montage grid (Figure 6): the **full 64-channel** montage is *confirmatory* at every
amplitude; the **ERP-CORE-30** montage falls to *exploratory* at 0.25°. The gate has
teeth precisely where the proposal feared a silent failure (sparse montage,
sub-degree movement).

### H2-analog — the headline audit, gated honestly (E5) ✅

A simulated ERP-CORE-like N2pc target (sparse montage, genuine posterior N2pc, plus a
**sub-0.2° condition-correlated ocular residual** of the kind that survives ERP CORE's
0.2° HEOG rejection): mean OPI **+0.150**, ceiling AUC 0.777, permutation p<0.005 —
**but the reliability gate returns "inconclusive"** (gaze R² 0.429 < t_expl). This is
the **Stage-0 "inconclusive by design" safeguard working**: a positive, significant
index is *withheld from certification* because the montage+regime cannot reconstruct
gaze well enough to trust it. The instrument refuses to over-claim against a trusted
dataset — exactly the behavior the proposal promises.

### Real-data negative control — PhysioNet eegmmidb (H3) ✅

Run on **real 64-channel EEG** (8 subjects downloaded). The ocular subspace is frozen
from 4 *source* subjects' eyes-open baseline and applied to 4 *held-out target*
subjects' motor runs (`real_eegmmidb_audit.csv`, Figure 7):

| contrast | ceiling AUC | ocular-proj AUC | OPI | p |
|---|---|---|---|---|
| **left vs right imagery** (no differential gaze) | 0.601 | 0.501 | **−0.023** | **0.998** |
| rest vs move (frontal-proxy) | 0.672 | 0.615 | +0.114 | 0.002 |

The **left-vs-right** contrast — where the eyes have no reason to differ — gives an
essentially perfect **null** on real EEG: the motor label is (weakly) decodable from
full EEG but **not** from the ocular subspace. The **rest-vs-move** contrast is
included deliberately as a cautionary result: a *frontal-channel-proxy* subspace
(the best one can build with **no eye-tracker in eegmmidb**) retains frontal *neural*
signal it cannot separate from ocular activity. That is not a bug — it is the
empirical argument for the proposal's central design choice: **the auditor must be
trained on synchronized gaze (EEGEyeNet), which isolates the ocular subspace from
frontal neural activity.** In the simulation, where gaze ground truth *is* available,
the analogous motor control collapses to OPI 0.030.

### Neural-leakage bound — leakage is governed by a measurable principal angle (E6) ✅

Added when the idea was strengthened. `R/run_leakage.R`
sweeps a lateralized neural source from posterior to frontal with **no ocular coupling**,
so any positive OPI is pure neural leakage, and relates it to the overlap
`s = ‖P t‖/‖t‖` between the neural contrast topography `t` and the ocular subspace `P`
(`results/tables/E6_leakage.csv`, Figure 8):

| | min overlap (posterior) | max overlap (frontal) | Spearman ρ(s, leakage) |
|---|---|---|---|
| overlap `s` | 0.02 | 0.6–0.8 | — |
| leakage OPI | **+0.006** | **+0.12 to +0.14** | **0.95** |

Leakage rises monotonically with `s` and is ≈0 when the cognitive topography is
orthogonal to the ocular subspace — which is exactly the posterior **N2pc** headline
case. This turns the auditor's central failure mode into a **per-dataset computable
number**: report `s` between the target's cognitive-decoder topography and the ocular
subspace, and a small `s` certifies low leakage. (Stated honestly: `s` upper-bounds the
overlap; the monotonicity is this pilot result under the linear-subspace model.)

## 4. Hypothesis scorecard

| | status in this run |
|---|---|
| **H0** Stage-0 feasibility gate ("inconclusive by design") | ✅ demonstrated (E4 gate + E5 headline gated inconclusive) |
| **H1** montage-matched reconstruction above threshold, both regimes | ✅ dense across regimes; ⚠️ sparse+fixational drops — *quantified, not assumed*. NB: real EEGEyeNet has **no** fixational regime (see FACTCHECK), so the near-fixation anchor must be ZuCo/N2pc, not EEGEyeNet |
| **H2** headline audit exceeds null by > calibrated margin | ✅ in-silico (E5); real ERP CORE audit is the compute-bound next step the shipped tool enables |
| **H3** graded sensitivity (pos high; neg null; covert tests leakage) | ✅ fully (E1, E2) + ✅ real data (eegmmidb L/R) |
| **H4** projection vs SP-overweighted ICLabel-gated ICA | ❌ not implemented — needs the OPTICAT/ICLabel Python/EEGLAB stack (future work) |

## 5. Honest scope — what is real, simulated, and not yet done

**Real:** the full auditor + locked projector + OPI + permutation null + reliability
gate (R, unit-tested); a real 64-ch EEG audit on PhysioNet eegmmidb (downloaded,
parsed, run) demonstrating the H3 negative control with a frozen cross-subject ocular
subspace.

**Simulated (with ground truth, by necessity):** the H3 battery, coupling sweep, bias
envelope, reliability calibration, and the H2-analog headline. Ground-truth gaze —
required to *measure* the bias envelope and prove the index recovers the true ocular
contribution — exists only in simulation (or in a synchronized-gaze dataset like
EEGEyeNet/ZuCo, not downloaded here for size/time).

**Not done (compute/data-bound, flagged honestly):** training the auditor on the full
356-subject EEGEyeNet corpus; the real ERP CORE N2pc/P3 audit; the ZuCo across-montage
bias anchor; and H4 (the ICA comparison). These need the large EEG datasets and the
Python EEG/DL stack (MNE/EEGLAB/OPTICAT); the shipped R tool is the instrument that
makes them one script away.

## 6. The shipped tool

`audit_dataset(X, montage, labels, subject, auditor)` in `R/auditor.R` is CPU-runnable
and returns the OPI, its ceiling, the permutation p, the two-sided envelope, and the
trust tier for any EEG matrix + montage — the deliverable the proposal promises
(Method 7).
