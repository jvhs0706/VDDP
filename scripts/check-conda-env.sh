#!/bin/bash
set -e

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