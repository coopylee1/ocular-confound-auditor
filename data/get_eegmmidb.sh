#!/usr/bin/env bash
# Download the PhysioNet eegmmidb subset used by R/run_eegmmidb.R.
# 8 subjects x 4 runs (R01 eyes-open baseline; R04/R08/R12 motor imagery) ~ 75 MB.
# Data are NOT redistributed in this repo; this fetches them from the original source.
set -e
cd "$(dirname "$0")"

BASE="https://physionet.org/files/eegmmidb/1.0.0"
SUBJECTS="S001 S002 S003 S004 S005 S006 S007 S008"
RUNS="R01 R04 R08 R12"

for s in $SUBJECTS; do
  mkdir -p "eegmmidb/$s"
  for r in $RUNS; do
    out="eegmmidb/$s/${s}${r}.edf"
    if [ -f "$out" ]; then
      echo "have   $out"
    else
      echo "get    $out"
      curl -fsSL "$BASE/$s/${s}${r}.edf" -o "$out"
    fi
  done
done

echo "Done. eegmmidb subset is in data/eegmmidb/"
