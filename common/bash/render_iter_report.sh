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
echo "SCRIPT_DIR: $SCRIPT_DIR"
echo "PROJECT_ROOT: $PROJECT_ROOT"
cd "$PROJECT_ROOT" || {
  echo "❌ ERROR: Could not cd to project root"
  exit 1
}

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
# ✅ Define output paths
# -------------------------------------
BASENAME="${EXP_ID}_${SIM_ID}_${ITER_ID}_report.html"
OUTPUT_DIR="docs/iter_reports"
FINAL_PATH="${OUTPUT_DIR}/${BASENAME}"

mkdir -p "$OUTPUT_DIR"

# -------------------------------
# ✅ Render to working directory
# -------------------------------
echo "🔧  Rendering report..."
quarto render "$QMD_PATH" \
  --to html-self-contained \
  --execute-params "$CONFIG_PATH" \
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
ln -s "$FINAL_PATH" "${ITER_DIR}/report.html"

echo "✅  Report saved to: $FINAL_PATH"
echo "🔗  Symlink created at: ${ITER_DIR}/report.html → $RELATIVE_LINK"
