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

NUM_REPEAT=5
COUNT_RANGE=114514

# if logs/vddlm.csv exists, remove it
if [ -f "logs/vddlm.csv" ]; then
    rm logs/vddlm.csv
fi

for i in $(seq 1 $NUM_REPEAT)
do 
    for dim in 4 16 64 256 1024 4096 16384 65536
    do
        for eps in 0.125 0.25 0.5 1.0 2.0 4.0 8.0
        do
            for noise_log_range in 4 6 8 10 12 14 16
            do
                for prec in 4 6 8 10 12 14 16
                do
                    python scripts/vddlm.py --dim $dim --eps $eps --count_range $COUNT_RANGE --noise_log_range $noise_log_range --prec $prec
                done
            done
        done
    done
done 