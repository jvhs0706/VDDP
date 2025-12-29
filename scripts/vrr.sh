#!/bin/bash
#SBATCH --time=23:59:59
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --job-name=VRR
#SBATCH --output=logs/vrr-%N-%j.out
#SBATCH --error=logs/vrr-%N-%j.err

#!/bin/bash
set -e
scripts/check-conda-env.sh

if [ ! -d "build" ]; then
    mkdir build
fi
cd build
cmake ..
make
cd ..

NUM_REPEAT=10

LOG_FILE="logs/vrr.csv"

# if logs/vrr.csv exists, remove it
if [ -f $LOG_FILE ]; then
    rm $LOG_FILE
fi

for i in $(seq 1 $NUM_REPEAT)
do 
    for eps in 0.125 0.25 0.5 1.0 2.0 4.0 8.0
    do
        for log_prec in 6 7 8 9 10 11 12
        do
            for num_class in 2 4 8 16 32 64 128
            do
                python scripts/vrr.py --eps $eps --log_prec $log_prec --num_class $num_class
            done
        done
        git add $LOG_FILE
        git commit -m "Update vrr logs for eps = $eps, run $i"
        git push
    done
done 