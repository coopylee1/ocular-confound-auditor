SUBMISSION PACKAGE
Title: Quantifying Ocular Confounds in EEG Decoding Without an Eye Tracker

Authors: Cooper Lee, Andrew Wong (corresponding author, aawong@umich.edu),
Ankit Walishetti, Luke Yin, Adam Park, and Yuan Yang. Affiliations appear on the
manuscript title page.

CONTENTS
  manuscript.tex          Main LaTeX source (IEEE conference template, IEEEtran;
                          self-contained, references inlined as \thebibliography)
  manuscript.pdf           Typeset manuscript (reference copy)
  fig4_coupling.png        Figure 1 (graded sensitivity vs. gaze-label coupling)
  fig8_leakage_angle.png   Figure 2 (neural-leakage diagnostic vs. subspace overlap)

COMPILATION
  pdflatex manuscript
  pdflatex manuscript

The bibliography is inlined in manuscript.tex (the IEEEtran.bst output, per the
IEEEtran HOWTO recommendation for submission), so no .bib/.bst/.bbl files and no
BibTeX run are required. Two pdflatex passes resolve all cross-references.

DATA AVAILABILITY
  All datasets analysed are publicly available from their original repositories
  (EEGEyeNet, ERP CORE, ZuCo, PhysioNet eegmmidb, Brain Invaders); see the
  manuscript's Acknowledgment section. No third-party data are redistributed in
  this package.
