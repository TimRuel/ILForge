#!/bin/bash
#SBATCH --account=p32397
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=timothyruel2024@u.northwestern.edu
#SBATCH --job-name=experiment_sim
#SBATCH --output=/dev/null
#SBATCH --nodes=1
#SBATCH --ntasks=1                 # one R process
#SBATCH --cpus-per-task=64         # all forked workers come from this
#SBATCH --mem=64G                  # total memory for the task
#SBATCH --array=0

# ===============================
# ✅ Validate CLI arguments
# ===============================
if [[ $# -ne 5 ]]; then
  echo "❌ ERROR: Missing arguments."
  echo "Usage: sbatch $0 <application_name> <estimand_name> <model_name> <experiment_id> <simulation_id>"
  exit 1
fi

APP_NAME="$1"
ESTIMAND="$2"
MODEL="$3"
EXP_ID="$4"
SIM_ID="$5"

ITER_NUM=$((SLURM_ARRAY_TASK_ID + 1))
ITER_ID=$(printf "iter_%04d" "$ITER_NUM")
REQUESTED_CORES=${SLURM_CPUS_PER_TASK:-1}

# ===============================
# ✅ Resolve project root
# ===============================
PROJECT_ROOT="$SLURM_SUBMIT_DIR"
cd "$PROJECT_ROOT" || {
  echo "❌ ERROR: Failed to cd into SLURM_SUBMIT_DIR ($SLURM_SUBMIT_DIR)"
  exit 1
}
echo "📁 PROJECT_ROOT resolved to: $PROJECT_ROOT"

# ===============================
# ✅ Logging setup
# ===============================
SIM_DIR="experiments/${EXP_ID}/simulations/${SIM_ID}"
ITER_DIR="${SIM_DIR}/${ITER_ID}"
LOG_DIR="${ITER_DIR}/logs"
mkdir -p "$LOG_DIR"

LOG_FILE="${LOG_DIR}/slurm_log.out"
CHECKJOB_MONITOR="${LOG_DIR}/checkjob_monitor.out"

exec > "$LOG_FILE" 2>&1
echo "📌 Logging to $LOG_FILE"

# ===============================
# ✅ Load environment modules
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
module load git-lfs/3.3.0-gcc-10.4.0

git lfs version || { echo "❌ git-lfs not available"; exit 1; }

# --- Limit BLAS threading conflicts (critical for multicore) ---
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "🔁 Running Iteration $ITER_NUM of Simulation $SIM_ID in Experiment $EXP_ID ($APP_NAME / $ESTIMAND/ $MODEL) with $REQUESTED_CORES cores..."

# ===============================
# ✅ Background checkjob monitor
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
# ✅ Cleanup diagnostics
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
# ✅ Run main R script
# ===============================
RSCRIPT_PATH="common/scripts/main.R"

if [[ ! -f "$RSCRIPT_PATH" ]]; then
  echo "❌ ERROR: Could not find R script at: $RSCRIPT_PATH"
  exit 1
fi

command -v Rscript >/dev/null 2>&1 || {
  echo "❌ ERROR: Rscript not found in PATH."
  exit 1
}

Rscript --max-connections=256 "$RSCRIPT_PATH" \
  "$APP_NAME" "$ESTIMAND" "$MODEL" "$EXP_ID" "$SIM_ID" "$ITER_ID" "$REQUESTED_CORES"

echo "✅ SLURM iteration complete: $ITER_ID"
