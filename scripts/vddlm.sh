#!/bin/bash
#SBATCH --time=6-23:59:59
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --job-name=VDDLM
#SBATCH --output=logs/vddlm-%N-%j.out
#SBATCH --error=logs/vddlm-%N-%j.err

#!/bin/bash
set -e
scripts/check-conda-env.sh
scripts/build.sh

DELTA_MIN=1e-10
LOG_NOISE_RANGE_SEARCH_SPACE=50
PREC_SEARCH_SPACE=50
COUNT_RANGE=114514

LOG_FILE="logs/vddlm.csv"

# if logs/vddlm.csv exists, remove it
if [ -f $LOG_FILE ]; then
    rm $LOG_FILE
fi

NUM_REPEAT=10

for i in $(seq 1 $NUM_REPEAT)
do
    for dim in 16 64 256 1024 4096
    do
        for eps in 1e-3 1e-2 1e-1 1e0
        do
            python scripts/vddlm.py --dim $dim --eps $eps --count_range $COUNT_RANGE --delta $DELTA_MIN --log_noise_range_search_space $LOG_NOISE_RANGE_SEARCH_SPACE --prec_search_space $PREC_SEARCH_SPACE 
            git add $LOG_FILE
            git commit -m "Update vddlm logs for dim $dim and eps $eps, run $i"
            git push
        done
    done 
done