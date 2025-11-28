#!/bin/bash
# This script exports the current conda environment to a YAML file.

# Check if conda is installed, and make sure the current environment is VDDP_ENV

if ! command -v conda &> /dev/null
then
    echo "Conda could not be found. Please install Conda and try again."
    exit 1
fi

CURRENT_ENV=$(conda info --json | jq -r '.active_prefix_name')
if [ "$CURRENT_ENV" != "VDDP_ENV" ]; then
    echo "Please activate the VDDP_ENV conda environment before running this script."
    exit 1
fi

# Export the current conda environment to a env.yaml, excluding the prefix line
conda env export --no-builds | grep -v "prefix: " > env.yaml
echo "Conda environment exported to env.yaml"