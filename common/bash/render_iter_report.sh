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
# ✅ Resolve project root
# -------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

# -------------------------------
# ✅ Construct config path
# -------------------------------
ITER_DIR="experiments/${EXP_ID}/simulations/${SIM_ID}/${ITER_ID}"
CONFIG_PATH="${ITER_DIR}/config_snapshot.yml"

if [ ! -f "$CONFIG_PATH" ]; then
  echo "❌ ERROR: Config file not found at: $CONFIG_PATH"
  exit 1
fi

# -------------------------------
# ✅ Find the .qmd template
# -------------------------------
QMD_PATH="applications/${APP}/${ESTIMAND}/quarto/iter_report.qmd"

if [ ! -f "$QMD_PATH" ]; then
  echo "❌ ERROR: Quarto template not found at: $QMD_PATH"
  exit 1
fi

# -------------------------------------
# ✅ Define output name and render
# -------------------------------------
OUTPUT_NAME="${EXP_ID}_${SIM_ID}_${ITER_ID}_report.html"
OUTPUT_DIR="docs/iter_reports"
OUTPUT_PATH="${OUTPUT_DIR}/${OUTPUT_NAME}"

mkdir -p "$OUTPUT_DIR"

quarto render "$QMD_PATH" \
  --execute-params "$CONFIG_PATH" \
  --output "$OUTPUT_NAME"

mv "$OUTPUT_NAME" "$OUTPUT_PATH"

# ----------------------------------------------
# ✅ Symlink report in iteration directory
# ----------------------------------------------
ln -sf "../../../../../$OUTPUT_PATH" "${ITER_DIR}/report.html"

echo "✅ Report saved to: $OUTPUT_PATH"
echo "🔗 Symlink created in: ${ITER_DIR}/report.html"
