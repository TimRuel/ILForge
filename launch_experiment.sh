#!/bin/bash
#SBATCH --account=p32397
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=timothyruel2024@u.northwestern.edu
#SBATCH --job-name=experiment_sim
#SBATCH --output=/dev/null
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=500M
#SBATCH --ntasks-per-node=64
#SBATCH --array=0-999

# ===============================
# ✅ Validate CLI Arguments
# ===============================
if [[ $# -lt 4 ]]; then
  echo "❌ ERROR: Missing arguments."
  echo "Usage: sbatch $0 <application_name> <estimand_name> <experiment_version> <sim_number>"
  exit 1
fi

APP_NAME="$1"               # e.g. multinomial_logistic_regression
ESTIMAND="$2"               # e.g. entropy
EXPERIMENT_VERSION="$3"     # e.g. v1.0.0
SIM_NUM="$4"                # e.g. 1

EXPERIMENT_ID=$(printf "exp_%s" "$EXPERIMENT_VERSION")
SIM_ID=$(printf "sim_%02d" "$SIM_NUM")
RUN_NUM=$((SLURM_ARRAY_TASK_ID + 1))
RUN_ID=$(printf "iter_%04d" "$RUN_NUM")
REQUESTED_CORES=$SLURM_NTASKS

SIM_DIR="experiments/${EXPERIMENT_ID}/simulations/${SIM_ID}"
RUN_DIR="${SIM_DIR}/${RUN_ID}"
LOG_DIR="${RUN_DIR}/logs"

mkdir -p "$LOG_DIR"

LOG_FILE="${LOG_DIR}/slurm_log.out"
CHECKJOB_MONITOR="${LOG_DIR}/checkjob_monitor.out"

# Redirect stdout and stderr to log file
exec > "$LOG_FILE" 2>&1

echo "📌 Logging to $LOG_FILE"

# ===============================
# ✅ Load Environment
# ===============================
module purge all
module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0 
module load gsl/2.7.1-gcc-12.3.0 
module load fftw/3.3.10-gcc-12.3.0 
module load gdal/3.7.0-gcc-12.3.0
module load nlopt/2.7.1-gcc-12.3.0
module load git/2.37.2-gcc-10.4.0
module load chrome/114.0.5735.90

echo "🔁 Running Iteration $RUN_NUM of Simulation $SIM_NUM in Experiment $EXPERIMENT_VERSION ($APP_NAME / $ESTIMAND) with $REQUESTED_CORES cores..."

# ===============================
# ✅ Monitor checkjob in background
# ===============================
(
  while true; do
    echo "===== checkjob (interval) at $(date) =====" >> "$CHECKJOB_MONITOR"
    checkjob "$SLURM_JOB_ID" >> "$CHECKJOB_MONITOR" 2>&1
    sleep 30
  done
) &
CHECKJOB_PID=$!

# ===============================
# ✅ Cleanup diagnostics on exit
# ===============================
trap '
  echo "===== FINAL DIAGNOSTICS for Job $SLURM_JOB_ID =====" >> "$LOG_FILE"
  seff "$SLURM_JOB_ID" >> "$LOG_FILE" 2>&1
  checkjob "$SLURM_JOB_ID" >> "$LOG_FILE" 2>&1
  sacct -j "$SLURM_JOB_ID" --format=JobID,JobName%20,Elapsed,MaxRSS,ReqMem,AllocCPUs,State >> "$LOG_FILE" 2>&1
  free -h >> "$LOG_FILE" 2>&1
  uptime >> "$LOG_FILE" 2>&1
  top -b -n 1 | head -40 >> "$LOG_FILE" 2>&1
  echo "===== BACKGROUND CHECKJOB MONITOR LOG =====" >> "$LOG_FILE"
  cat "$CHECKJOB_MONITOR" >> "$LOG_FILE" 2>/dev/null
  rm -f "$CHECKJOB_MONITOR"
  kill "$CHECKJOB_PID" 2>/dev/null
' EXIT

# ===============================
# ✅ Run R script
# ===============================
RSCRIPT_PATH="applications/${APP_NAME}/${ESTIMAND}/scripts/main.R"

if [[ ! -f "$RSCRIPT_PATH" ]]; then
  echo "❌ ERROR: Could not find R script at: $RSCRIPT_PATH"
  exit 1
fi

command -v Rscript >/dev/null 2>&1 || {
  echo "❌ ERROR: Rscript not found in PATH."
  exit 1
}

EXTRA_FLAGS=("--force" "--non-interactive")

Rscript --max-connections=256 "$RSCRIPT_PATH" \
  "$APP_NAME" "$ESTIMAND" "$EXPERIMENT_ID" "$REQUESTED_CORES" "$SIM_ID" "$RUN_ID" \
  --force --non-interactive # Add: --skip-integrated, --skip-profile if needed


