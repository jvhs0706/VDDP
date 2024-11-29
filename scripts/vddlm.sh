#!/bin/bash
#SBATCH --time=2-23:59:59
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

if [ ! -d "build" ]; then
    mkdir build
fi
cd build
cmake ..
make
cd ..

NUM_REPEAT=1
DELTA_MIN=1e-10
LOG_NOISE_RANGE_SEARCH_SPACE=50
PREC_SEARCH_SPACE=50
COUNT_RANGE=114514

# if logs/vddlm.csv exists, remove it
if [ -f "logs/vddlm.csv" ]; then
    rm logs/vddlm.csv
fi

for dim in 16 64 256 1024 4096
do
    for eps in 9.5e-4 9.5e-3 9.5e-2 9.5e-1
    do
        for i in $(seq 1 $NUM_REPEAT)
        do
            python scripts/vddlm.py --dim $dim --eps $eps --count_range $COUNT_RANGE --delta $DELTA_MIN --log_noise_range_search_space $LOG_NOISE_RANGE_SEARCH_SPACE --prec_search_space $PREC_SEARCH_SPACE 
        done
    done
done 