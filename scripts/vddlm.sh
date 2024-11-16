source activate ~/.conda/envs/VDP_ENV
echo "Using cmake $(which cmake), version $(cmake --version)."

if [ ! -d "build" ]; then
    mkdir build
fi
cd build
cmake ..
make
cd ..

./build/src/vddlm 1024 0.5 6 10