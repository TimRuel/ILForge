#!/usr/bin/bash

# Usage: ./common/bash/render_iter_report.sh <APP> <ESTIMAND> <EXP_ID> <SIM_ID> <ITER_ID>
# Example: ./common/bash/render_iter_report.sh multinomial_logistic_regression entropy exp_v1.0.0 sim_01 iter_0001

set -e

# -------------------------------
# ✅ Validate arguments
# -------------------------------
if [ "$#" -ne 5 ]; then
  echo "Usage: $0 <APP> <ESTIMAND> <EXP_ID> <SIM_ID> <ITER_ID>"
  exit 1
fi

APP="$1"
ESTIMAND="$2"
EXP_ID="$3"
SIM_ID="$4"
ITER_ID="$5"

# -------------------------------
# ✅ Resolve project root and QMD directory
# -------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
QMD_DIR="${PROJECT_ROOT}/applications/${APP}/${ESTIMAND}/quarto"
QMD_BASENAME="iter_report.qmd"
QMD_PATH="${QMD_DIR}/${QMD_BASENAME}"

# -------------------------------
# ✅ Validate input files
# -------------------------------
ITER_DIR="${PROJECT_ROOT}/experiments/${EXP_ID}/simulations/${SIM_ID}/${ITER_ID}"
CONFIG_PATH="${ITER_DIR}/config_snapshot.yml"

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
BASENAME="${EXP_ID}_${SIM_ID}_${ITER_ID}_report.html"
OUTPUT_DIR="${PROJECT_ROOT}/docs/iter_reports"
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
# ✅ Create symlink in iteration folder
# -------------------------------
RELATIVE_LINK="../../../../../docs/iter_reports/${BASENAME}"
ln -sf "$RELATIVE_LINK" "${ITER_DIR}/iter_report.html"

echo "✅ Report saved to: $FINAL_PATH"
echo "🔗 Symlink created at: ${ITER_DIR}/iter_report.html → $FINAL_PATH"
