# Camera-ready changes — paper #1571313658, IEEE Healthcom 2026

*Quantifying Ocular Confounds in EEG Decoding Without an Eye Tracker*

Every number below was produced by running code in `R/`; the new analyses are in
`R/run_revision.R`, `R/run_revision2.R`, `R/fig_revision.R`, with outputs in
`results/tables/E7*.csv`, `E8*.csv`, `E9*.csv`, `E10*`.

---

## The one thing we did NOT do

**Reviewers 1, 2 and 3 all asked, as their top request, for the auditor to be trained on real
synchronised EEG–gaze data (EEGEyeNet / ZuCo) and tested on independent participants. That is not
in this camera-ready.** It is a new-data experiment, not a revision, and we judged it
undeliverable at camera-ready quality in the time available. (The EEGEyeNet OSF record's Dropbox
provider is dead; the data now sits behind a Google Drive provider and is multi-GB.)

Instead of claiming it, we have **narrowed the paper's claims to match the evidence** and stated
the missing experiment explicitly as the decisive next step, in the abstract, the limitations,
and the conclusion. The paper is now presented as a *simulation-validated proof of concept with a
preliminary real-EEG illustration*, which is precisely R1's recommended framing.

---

## Terminology and claim scope

| Reviewer request | Change |
|---|---|
| R1: "ocular subspace" → "gaze-predictive EEG subspace" | Renamed throughout, with a Methods paragraph explaining *why* the rename is required (the readout may capture gaze-correlated cortical activity). |
| R1, R3: OPI is not a fraction or causal estimate | Redefined in the abstract and in a dedicated Methods paragraph as *label predictability within the projection*. Explicitly states that a high OPI does not establish ocular causation and a low OPI does not clear a dataset (may reflect failed transfer). Removed all "how much of a result is ocular" phrasing. |
| R1: avoid "bound" for the overlap score | Replaced with "diagnostic" in text **and in the figure's own embedded title**, which still said "bounded" in the submitted version. |
| R1: field-wide deployment is premature | Deleted the leaderboard proposal. The paper now says explicitly that it proposes no cross-benchmark leaderboard on this evidence. |
| R1: narrow to simulation-validated proof of concept | Done in abstract, discussion, limitations, conclusion. |

## New analyses

**Baselines (R1 ask #10) — Table I.** OPI computed under four projections on identical data:
the gaze-trained subspace (k=4), a **rank-matched random** subspace (k=4), a bipolar **EOG proxy**
(k=2, F8−F7 and (Fp1+Fp2)/2), and a **frontal-channel audit** (k=17).

This changed a conclusion. At the 8° saccade amplitude used for the original headline battery,
*every* projection — including the random one — returns OPI ≈ 0.43. The headline positive control
therefore carries no information about the merit of gaze training; at 8° the artifact simply
dominates everything. The methods only separate in the fixational regime:

| amplitude | gaze | random | EOG | frontal |
|---|---|---|---|---|
| 0.25° | **0.396** | 0.149 | 0.144 | 0.401 |
| 8° | 0.429 | 0.428 | 0.429 | 0.429 |

On the covert-attention neural control the random projection returns OPI 0.072 and is significant
in **92%** of studies — a false ocular attribution — where the gaze-trained subspace returns 0.007
(17%). We report honestly that the **frontal-channel audit matches the gaze-trained operator on
both axes in this simulator**, and that we therefore claim no superiority over it; the simulator
contains no frontal neural generator, which is exactly what the real-data proxy experiment shows
it would need.

**Per-output reliability (R1 ask #8) — Table II.** R² reported separately for horizontal gaze,
vertical gaze and the saccade channel; **AUC** for the binary blink indicator (an R² on a 0/1
outcome is not interpretable). This also changed a conclusion: the aggregate was propped up by
near-perfect blink detection and easy vertical gaze, hiding that *horizontal* gaze — the
dimension a lateralized paradigm depends on — collapses at the fixational extreme (R² 0.474 full,
0.353 on 30 channels). Under the aggregate the full montage at 0.25° passes the confirmatory
threshold (0.664); under the horizontal-gaze score it does not (0.474). The per-output gate is
strictly more conservative and reverses a verdict.

**Held-out threshold calibration (R1 ask #7) — §III-D.** Thresholds were previously fit and
evaluated on the same grid. Refit on the full-64 grid alone and applied unchanged to the held-out
30-channel grid: confirmatory threshold R² ≈ 0.66, six of seven held-out cells tiered
confirmatory, all six with |bias| ≤ 0.05, and the 0.25° cell correctly excluded. Under honest
calibration **the confirmatory and exploratory thresholds coincide, so the evidence supports a
two-tier gate, not the three tiers previously reported.** Values rounded to 2 s.f. and labelled
preliminary.

**Activation pattern vs weight vector (R1 ask #9) — §III-E.** We computed the Haufe transform
a = Σ_X t as requested. It is *counter-productive here*, and we report why: Σ_X is dominated by
the ocular variance the diagnostic is meant to be separate from, so a points into the
gaze-predictive subspace almost regardless of t. Its overlap saturates at 0.92–0.99 across every
source position and correlates **negatively** with leakage (ρ = −0.54), where the discriminative
direction gives ρ = +0.95 and ranks conditions *identically* to the true source topography
(ρ = 1.00). The paper now argues explicitly that the activation pattern is the right object for
asking what generated a component and the wrong object for asking what a decoder reads.

**Uncertainty (R1 ask #11, R3).** 95% CIs on all headline OPI values, CI bars on the figure, SDs
noted for the reliability table. Thresholds no longer quoted to three decimals.

## Motivation (R2)

R2 could not tell whether this addresses a real problem or is "a nice mathematical exercise." A
new Introduction paragraph gives the concrete clinical stake: a communication BCI whose decodable
signal is substantially ocular fails precisely in ALS and other oculomotor-impaired populations —
the patients it targets; a gaze-driven attention/workload biomarker tracks looking behaviour and
will not transfer across populations with different eye-movement statistics. The abstract now
opens on this.

R2 also asked for evaluation against a real eye-tracker. That is the experiment we did not do;
see the top of this document.

## Real-data framing (R1, R3)

The PhysioNet analysis is now labelled an **illustrative proxy** in its section title, methods,
and results, with explicit statements that (a) it does not implement the proposed gaze-trained
operator, (b) the left/right no-differential-gaze assumption is *assumed, not measured*, since the
dataset has no eye-tracker, and (c) eight participants split four/four is too small for stable
estimates. Its value is stated as negative evidence.

## Trims for the 6-page limit

Condensed physiology background, removed the repeated full-channel vs 30-channel discussion,
deleted the field-wide leaderboard paragraph, moved simulator parameters to the released code, and
cut 5 references. **Two items were cut purely for space and can be restored if Healthcom permits a
7th page:** the gaze–label coupling figure (its six values are now given inline in the text) and
the real-EEG results table (now given inline in prose).

## Files

- `manuscript.tex` / `manuscript.pdf` — camera-ready, 6 pages
- `fig8_leakage_angle.png` — regenerated with corrected terminology and CI bars
- New analysis code: `R/run_revision.R`, `R/run_revision2.R`, `R/fig_revision.R`
