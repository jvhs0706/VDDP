#!/bin/bash
scripts/check-conda-env.sh

if [ ! -d "build" ]; then
    mkdir build
fi
cd build
cmake ..
make
cd ..