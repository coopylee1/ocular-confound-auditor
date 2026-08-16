CAMERA-READY SUBMISSION PACKAGE
IEEE Healthcom 2026 -- paper #1571313658

Title: Quantifying Ocular Confounds in EEG Decoding Without an Eye Tracker

Authors: Cooper Lee, Andrew Wong (corresponding author, aawong@umich.edu),
Ankit Walishetti, Luke Yin, Adam Park, and Yuan Yang. Affiliations appear on the
manuscript title page.

CONTENTS
  manuscript.tex           Main LaTeX source (IEEE conference template, IEEEtran;
                           self-contained, references inlined as \thebibliography)
  manuscript.pdf           Typeset manuscript (6 pages)
  fig8_leakage_angle.png   Figure 1 (neural leakage vs. subspace overlap)
  RESPONSE_TO_REVIEWERS.md Summary of camera-ready changes (not for publication)

COMPILATION
  pdflatex manuscript
  pdflatex manuscript

The bibliography is inlined in manuscript.tex (the IEEEtran.bst output, per the
IEEEtran HOWTO recommendation for submission), so no .bib/.bst/.bbl files and no
BibTeX run are required. Two pdflatex passes resolve all cross-references.

BEFORE UPLOADING -- author checklist
  [ ] At least one author registered and paid (registration MUST precede upload)
  [ ] PDF validated through IEEE PDF eXpress
  [ ] IEEE electronic copyright form (eCF) completed in EDAS
  [ ] Confirm the 6-page limit and whether over-length pages are available; two
      items were cut purely for space and can be restored if a 7th page is
      permitted (see RESPONSE_TO_REVIEWERS.md, "Trims for the 6-page limit")

DATA AVAILABILITY
  All datasets analysed are publicly available from their original repositories
  (EEGEyeNet, ERP CORE, ZuCo, PhysioNet eegmmidb); see the manuscript's
  Acknowledgment section. No third-party data are redistributed in this package.
