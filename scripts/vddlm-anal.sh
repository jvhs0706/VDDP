#!/bin/bash
#SBATCH --time=0:59:59
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=236511M
#SBATCH --job-name=AbwehrFuckAwayIfYouDontRespectMe
#SBATCH --output=logs/vddlm-%N-%j.out
#SBATCH --error=logs/vddlm-%N-%j.err
#SBATCH --exclude=snorlax-[1,3]

#SBATCH --mail-user=haochen.sun@uwaterloo.ca
#SBATCH --mail-type=ALL

source activate ~/.conda/envs/VDP_ENV

if [ -f "logs/vddlm-anal.csv" ]; then
    rm logs/vddlm-anal.csv
fi

DELTA_MIN=1e-10
LOG_NOISE_RANGE_SEARCH_SPACE=50
PREC_SEARCH_SPACE=50

for eps in 0.0009765625 0.001953125 0.00390625 0.0078125 0.015625 0.03125 0.0625 0.125 0.25 0.5 1.0
do
    python scripts/vddlm.py --eps $eps --delta $DELTA_MIN --log_noise_range_search_space $LOG_NOISE_RANGE_SEARCH_SPACE --prec_search_space $PREC_SEARCH_SPACE --anal
done

python scripts/vddlm-anal-plot.py