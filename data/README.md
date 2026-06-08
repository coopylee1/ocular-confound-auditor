# Data

This repository does **not** redistribute any third-party EEG data.

The real-data negative control (`R/run_eegmmidb.R`) uses the **EEG Motor Movement/Imagery
Dataset (`eegmmidb`)** from PhysioNet — 64-channel BCI2000 EEG, 160 Hz, EDF+:

> Schalk, G., McFarland, D.J., Hinterberger, T., Birbaumer, N., Wolpaw, J.R. (2004).
> BCI2000: A General-Purpose Brain-Computer Interface (BCI) System. *IEEE TBME* 51(6).
> Hosted on PhysioNet (Goldberger et al., 2000): https://physionet.org/content/eegmmidb/

It is openly available under the PhysioNet/ODC-BY terms. Download the exact subset this
project uses (8 subjects × 4 runs ≈ 75 MB) with:

```bash
bash data/get_eegmmidb.sh
```

This populates `data/eegmmidb/S001/…S008/` with runs `R01` (eyes-open baseline used to
estimate the frozen ocular subspace), `R04`, `R08`, `R12` (motor-imagery runs that are
audited). After downloading, `R/run_eegmmidb.R` (and `run_all.sh` step 5) will run.
