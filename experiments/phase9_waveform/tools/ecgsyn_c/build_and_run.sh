#!/bin/bash
# Rebuild the authors' ECGSYN C program (PhysioNet ecgsyn 1.0.0, GPL; not committed here) with the
# two stand-ins in this directory for the Numerical Recipes files dfour1.c / ran1.c (which only
# feed the RR process and the optional additive noise), and produce the reference runs used by
# experiments/phase9_waveform/verify_ecgsyn.py.
set -e
OUT=${1:-/tmp/claude-0/ecgsyn_build}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT" && cd "$OUT"
curl -sS -L "https://physionet.org/files/ecgsyn/1.0.0/ecgsyn.tar.gz?download" -o ecgsyn.tar.gz
tar xzf ecgsyn.tar.gz C/src/ecgsyn.c C/src/opt.c C/src/opt.h
cp C/src/ecgsyn.c C/src/opt.c C/src/opt.h "$HERE/dfour1.c" "$HERE/ran1.c" .
gcc -O -o ecgsyn ecgsyn.c opt.c dfour1.c ran1.c -lm
for spec in "60 0 r60c 60" "75 0 r75 120" "60 5 r60v 120" "90 3 r90v 120"; do
  set -- $spec; mkdir -p $3; (cd $3 && ../ecgsyn -n $4 -s 360 -S 720 -h $1 -H $2 -a 0 -R 7 -O ecg.dat </dev/null > log.txt)
done
echo "reference runs in $OUT"
