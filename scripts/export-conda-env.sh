#!/bin/bash
set -e
./check-conda-env.sh

# Export the current conda environment to a env.yaml, excluding the prefix line
conda env export --no-builds | grep -v "prefix: " > env.yaml
echo "Conda environment exported to env.yaml"