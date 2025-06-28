#!/bin/bash

# Usage: ./render_simulation_report.sh path/to/config_snapshot.yml
# Example: ./render_simulation_report.sh experiments/exp_v1.0.0/simulations/sim_01/iter_001/config_snapshot.yml

set -e

# Check for input argument
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 path/to/config_snapshot.yml"
  exit 1
fi

INPUT_PATH="$1"

# Resolve absolute path to config_snapshot.yml
if command -v realpath >/dev/null 2>&1; then
  CONFIG_PATH=$(realpath "$INPUT_PATH")
else
  CONFIG_PATH=$(cd "$(dirname "$INPUT_PATH")" && pwd)/$(basename "$INPUT_PATH")
fi

# Extract the parent directory of the directory containing config_snapshot.yml
ITER_DIR=$(dirname "$CONFIG_PATH")
SIM_DIR=$(dirname "$ITER_DIR")

# Define output filename
OUTPUT_FILE="simulation_report.html"

# Render the report
quarto render simulation_report.qmd \
  --execute-params "$CONFIG_PATH" \
  --output "$OUTPUT_FILE"

# Move the rendered file to the simulation directory
mv "$OUTPUT_FILE" "$SIM_DIR/"

echo "Simulation report saved to: $SIM_DIR/$OUTPUT_FILE"
