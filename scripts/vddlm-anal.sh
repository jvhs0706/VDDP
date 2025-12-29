#!/bin/bash
#SBATCH --time=0:59:59
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --job-name=VDDLM-ANAL
#SBATCH --output=logs/vddlm-%N-%j.out
#SBATCH --error=logs/vddlm-%N-%j.err

#!/bin/bash
set -e
./check-conda-env.sh

LOG_FILE="logs/vddlm-anal.csv"

if [ -f $LOG_FILE ]; then
    rm $LOG_FILE
fi

DELTA_MIN=1e-10
LOG_NOISE_RANGE_SEARCH_SPACE=50
PREC_SEARCH_SPACE=50

for eps in 0.0009765625 0.001953125 0.00390625 0.0078125 0.015625 0.03125 0.0625 0.125 0.25 0.5 1.0
do
    python scripts/vddlm.py --eps $eps --delta $DELTA_MIN --log_noise_range_search_space $LOG_NOISE_RANGE_SEARCH_SPACE --prec_search_space $PREC_SEARCH_SPACE --anal
done

python scripts/vddlm-anal-plot.py vddlm-anal

git add $LOG_FILE plots/vddlm-anal-plot.pdf
git commit -m "Update vddlm-anal logs"
git push