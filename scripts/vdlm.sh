source activate ~/.conda/envs/VDP_ENV

cd build
cmake ..
make
cd ..

OMP_NUM_THREADS=32 ./build/src/vdlm 1000 0.9 16