#!/bin/bash
#SBATCH --time=11:59:59
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=236511M
#SBATCH --job-name=shilishanlu
#SBATCH --output=logs/vrr-%N-%j.out
#SBATCH --error=logs/vrr-%N-%j.err

#SBATCH --mail-user=haochen.sun@uwaterloo.ca
#SBATCH --mail-type=ALL

source activate ~/.conda/envs/VDP_ENV

cd build
cmake ..
make
cd ..

NUM_REPEAT=100

# if logs/vrr.csv exists, remove it
if [ -f "logs/vrr.csv" ]; then
    rm logs/vrr.csv
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
    done
    python scripts/vrrplot.py # plot the results up to the current repeat
done 

