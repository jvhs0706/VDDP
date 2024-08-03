source activate ~/.conda/envs/VDP_ENV

cd build
cmake ..
make
cd ..

./build/src/vdlm 65536 1.0 4 16