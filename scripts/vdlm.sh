source activate ~/.conda/envs/VDP_ENV

cd build
cmake ..
make
cd ..

./build/src/vdlm 256 1.0 4 16