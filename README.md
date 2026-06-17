# Quantifying Ocular Confounds in EEG Decoding Without an Eye Tracker

A working **R** implementation and validation of the **Ocular Predictability Index
(OPI)** auditor: a frozen linear EEG→gaze readout, redeployed as an orthogonal
**ocular-subspace projection**, that measures how much of a cognitive-decoding result
*could be the eyes* — and that returns **"inconclusive"** rather than a false clean bill
of health when it cannot reconstruct gaze in the target montage.

The auditor operates on **all available channels** (a full 64-channel montage as primary,
with the 30-channel ERP-CORE set as an explicitly specified subset), is calibrated against
ground truth in a physically grounded EEG+gaze simulator, and is demonstrated on real
64-channel EEG (PhysioNet `eegmmidb`). See [`REPORT.md`](REPORT.md) for the full execution
report and [`paper/manuscript.pdf`](paper/manuscript.pdf) for the manuscript.

## Why R (and why that is faithful)
The **locked operator** is the orthogonal projection onto the span of the frozen **linear**
gaze/ocular-state readout weights, and the index is AUC-of-a-decoder on that projection
minus a within-subject permutation null. That is linear algebra + statistics — the object
that gets *frozen and shipped* is linear — so R implements the *actual locked instrument*,
not a stand-in. (A deep EEGEyeNet network is only one way to *train* the readout; the
shipped projector is linear either way.)

## Quickstart
```bash
git clone https://github.com/coopylee1/ocular-confound-auditor.git
cd ocular-confound-auditor
bash run_all.sh            # unit tests, demo, simulation, leakage, real audit, figures (~5 min, CPU)
Rscript R/demo_tool.R      # minimal usage of the shipped auditor
```
Dependencies (auto-installed once): `pROC, data.table, edfReader, signal, jsonlite, ggplot2`.

The simulation, leakage, and figure steps need no external data. The **real-data audit**
(`R/run_eegmmidb.R`) needs the PhysioNet `eegmmidb` EDF+ files under `data/eegmmidb/`,
which are **not redistributed here** — fetch them with:
```bash
bash data/get_eegmmidb.sh   # downloads the 8 subjects / 4 runs used (~75 MB) from PhysioNet
```
See [`data/README.md`](data/README.md).

## Layout
```
R/auditor.R         core library: ridge readout, LOCKED ocular projection,
                    OPI + within-subject permutation null, reliability gate,
                    audit_dataset() (the shipped tool)
R/simulate.R        grounded EEG+gaze forward model (corneo-retinal HEOG ~16uV/deg,
                    Keren saccadic spike, blink, N2pc); full 64-ch montage + ERP-CORE-30 subset
R/test_auditor.R    unit tests (AUC vs pROC, projector algebra, ridge, end-to-end)
R/run_simulation.R  E1 H3 battery · E2 coupling sweep · E3 bias envelope ·
                    E4 reliability calibration · E5 H2-analog headline
R/run_leakage.R     E6: neural-leakage term vs principal-angle overlap
R/run_eegmmidb.R    REAL 64-ch EEG negative control (frozen cross-subject ocular subspace)
R/figures.R         figures 1–9
R/demo_tool.R       minimal end-to-end usage of the shipped auditor
results/tables/*.csv, results/figures/*.png   all outputs (committed)
paper/              manuscript (.tex/.pdf), bibliography, figures
```

## Headline results (all produced by running the code; see REPORT.md)
- **Graded sensitivity (H3):** OPI **0.43** on saccade-coupled ocular contamination and
  **0.20** on a distinct condition-correlated **eye-blink** artifact (both 100% significant),
  versus **≤0.03** on every negative control (covert attention, motor imagery, label shuffle,
  neutral); rises monotonically with gaze↔label coupling (0 → 0.5).
- **Neural-leakage specificity:** a lateralized source decodable from full EEG (AUC **0.79**)
  is **not** decodable from the ocular projection (OPI **0.007**, n.s.) — the auditor does not
  absorb lateralized posterior neural signal.
- **Bias envelope:** auditor OPI matches the ground-truth-gaze OPI within **±0.01**, except a
  **conservative** −0.034 (full 64-ch) / −0.056 (30-ch) at the fixational extreme (0.25°) —
  measured, not assumed.
- **Reliability gate:** thresholds derived empirically (t_conf = **0.664**, t_expl = **0.497**);
  the 30-ch montage drops to *exploratory* at 0.25°, and the H2-analog headline is correctly
  gated **"inconclusive"** rather than over-claiming.
- **Neural-leakage bound:** leakage is governed by the principal-angle overlap between the
  cognitive topography and the ocular subspace (Spearman **ρ = 0.95**; ≈0 for a posterior source).
- **Real EEG (`eegmmidb`):** left-vs-right motor imagery → **OPI −0.023, p = 0.998** (clean null
  on real data). Rest-vs-move is included as a cautionary frontal-proxy-leakage case (+0.114).

## Scope (honest)
**Real:** the auditor + the real `eegmmidb` audit. **Simulated** (where ground-truth gaze is
*required* to measure bias): the H3 battery, reliability calibration, bias envelope, and the
H2-analog. **Not done** (data/compute-bound, flagged in `REPORT.md` §5): full 356-subject
EEGEyeNet training, the real ERP CORE audit, the ZuCo bias anchor, and the ICA comparison.

## Data availability
All datasets are publicly available from their original repositories (EEGEyeNet, ERP CORE,
ZuCo, PhysioNet `eegmmidb`, Brain Invaders). No third-party data are redistributed in this
repository; `data/get_eegmmidb.sh` fetches the `eegmmidb` runs used directly from PhysioNet.
