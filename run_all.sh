#!/usr/bin/env bash
# Reproduce the entire Paper-15 instrument + validation. CPU-only, ~5 min.
set -e
cd "$(dirname "$0")"
echo "== 1/5 unit tests ==";            Rscript R/test_auditor.R
echo "== 2/5 shipped-tool demo ==";     Rscript R/demo_tool.R
echo "== 3/6 in-silico validation ==";  Rscript R/run_simulation.R
echo "== 4/6 neural-leakage bound ==";  Rscript R/run_leakage.R
echo "== 5/6 real eegmmidb audit ==";   Rscript R/run_eegmmidb.R || echo "(skipped: download data/eegmmidb first)"
echo "== 6/6 figures ==";               Rscript R/figures.R
echo "Done. See results/tables and results/figures."
