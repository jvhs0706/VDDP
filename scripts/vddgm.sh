#!/bin/bash
#SBATCH --time=6-23:59:59
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --job-name=VDDGM
#SBATCH --output=logs/vddgm-%N-%j.out
#SBATCH --error=logs/vddgm-%N-%j.err

#!/bin/bash
set -e
./check-conda-env.sh

if [ ! -d "build" ]; then
    mkdir build
fi
cd build
cmake ..
make
cd ..

LOG_FILE="logs/vddgm.csv"

# if logs/vddgm.csv exists, remove it
if [ -f $LOG_FILE ]; then
    rm $LOG_FILE
fi

NUM_REPEAT=5

for i in $(seq 1 $NUM_REPEAT)
do
    for dim in 16 64 256 1024
    do
        for eps in 1 4 16 64
        do 
            python scripts/vddgm.py --dim $dim --sigma $eps
            git add $LOG_FILE
            git commit -m "Update vddgm logs for dim $dim and sigma $eps, run $i"
            git push
        done
    done
done