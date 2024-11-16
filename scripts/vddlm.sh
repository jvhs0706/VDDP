source activate ~/.conda/envs/VDP_ENV
echo "Using cmake $(which cmake), version $(cmake --version)."

if [ ! -d "build" ]; then
    mkdir build
fi
cd build
cmake ..
make
cd ..

python scripts/vddlm.py --dim 1024 --eps 1 --count_range 100 --noise_log_range 6 --prec 4