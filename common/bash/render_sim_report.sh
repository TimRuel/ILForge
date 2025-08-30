#!/usr/bin/bash

# Usage: ./common/bash/render_sim_report.sh <APP> <ESTIMAND> <MODEL> <EXP_ID> <SIM_ID>
# Example: ./common/bash/render_sim_report.sh poisson group_rates_weighted_sum naive_group_rates exp_v1.0.0 sim_01

set -e

# -------------------------------
# ✅ Validate arguments
# -------------------------------
if [ "$#" -ne 5 ]; then
  echo "Usage: $0 <APP> <ESTIMAND> <MODEL> <EXP_ID> <SIM_ID>" 
  exit 1
fi

APP="$1"
ESTIMAND="$2"
MODEL="$3"
EXP_ID="$4"
SIM_ID="$5"

# -------------------------------
# ✅ Resolve project root and QMD directory
# -------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
QMD_DIR="${PROJECT_ROOT}/applications/${APP}/${ESTIMAND}/${MODEL}/quarto"
QMD_BASENAME="sim_report.qmd"
QMD_PATH="${QMD_DIR}/${QMD_BASENAME}"

# -------------------------------
# ✅ Validate input files
# -------------------------------
SIM_DIR="${PROJECT_ROOT}/experiments/${EXP_ID}/simulations/${SIM_ID}"
CONFIG_PATH="${PROJECT_ROOT}/config/exps/${EXP_ID}.yml"

if [ ! -f "$CONFIG_PATH" ]; then
  echo "❌ ERROR: Config file not found at: $CONFIG_PATH"
  exit 1
fi

if [ ! -f "$QMD_PATH" ]; then
  echo "❌ ERROR: Quarto template not found at: $QMD_PATH"
  exit 1
fi

# -------------------------------
# ✅ Define output paths
# -------------------------------
BASENAME="${EXP_ID}_${SIM_ID}.html"
OUTPUT_DIR="${PROJECT_ROOT}/docs/sim_reports"
FINAL_PATH="${OUTPUT_DIR}/${BASENAME}"

mkdir -p "$OUTPUT_DIR"

# -------------------------------
# ✅ Render self-contained in QMD directory
# -------------------------------
echo "🔧 Rendering report..."
cd "$QMD_DIR"

# Make config path relative to QMD directory
CONFIG_RELATIVE="$(realpath --relative-to="$QMD_DIR" "$CONFIG_PATH")"

quarto render "$QMD_BASENAME" \
  -P sim_id:"$SIM_ID" \
  --execute-params "$CONFIG_RELATIVE" \
  --output "$BASENAME"

# -------------------------------
# ✅ Move to final location
# -------------------------------
if [ ! -f "$BASENAME" ]; then
  echo "❌ ERROR: Expected report not found: $BASENAME"
  exit 1
fi

mv "$BASENAME" "$FINAL_PATH"

# -------------------------------
# ✅ Create symlink in simulation folder
# -------------------------------
RELATIVE_LINK="../../../../docs/sim_reports/${BASENAME}"
ln -sf "$RELATIVE_LINK" "${SIM_DIR}/sim_report.html"

echo "✅ Report saved to: $FINAL_PATH"
echo "🔗 Symlink created at: ${SIM_DIR}/sim_report.html → $FINAL_PATH"
