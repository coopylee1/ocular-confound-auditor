SUBMISSION PACKAGE
Title: Quantifying Ocular Confounds in EEG Decoding Without an Eye Tracker

Authors: Cooper Lee, Andrew Wong (corresponding author, aawong@umich.edu),
Ankit Walishetti, Luke Yin, Adam Park, and Yuan Yang. Affiliations appear on the
manuscript title page.

CONTENTS
  manuscript.tex          Main LaTeX source (IEEE conference template, IEEEtran;
                          self-contained, references inlined as \thebibliography)
  manuscript.pdf           Typeset manuscript (reference copy)
  fig9_efield.png          Figure 1 (electric field of an eye movement)
  fig4_coupling.png        Figure 2 (graded sensitivity vs. gaze-label coupling)
  fig8_leakage_angle.png   Figure 3 (neural-leakage bound vs. subspace overlap)

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
