source activate ~/.conda/envs/VDP_ENV

cd build
cmake ..
make
cd ..

./build/src/vdlm 1048576 0.9 16