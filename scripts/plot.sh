#!/bin/bash
set -e
scripts/check-conda-env.sh

python scripts/vddlm-anal-plot.py vddlm-anal
python scripts/vddlmplot.py vddlm external-vdbm
python scripts/vddgmplot.py vddgm
python scripts/vrrcompplot.py vrr external-kcy21
python scripts/vrrdecompplot.py vrr

git add plots/*.pdf
git commit -m "Update plots"
git push