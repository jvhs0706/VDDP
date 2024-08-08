source activate ~/.conda/envs/VDP_ENV

cd build
cmake ..
make
cd ..

./build/src/vdlm 4096 1.0 4 16