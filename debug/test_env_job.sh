#!/bin/bash
#SBATCH --account=p32397
#SBATCH --partition=short
#SBATCH --time=00:05:00
#SBATCH --job-name=test_env
#SBATCH --output=test_env_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=500M

echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo "SLURM_JOB_NODELIST: $SLURM_JOB_NODELIST"
echo "SLURM_NTASKS: $SLURM_NTASKS"
echo "Working directory: $SLURM_SUBMIT_DIR"

# Load modules
module purge all
module load R/4.4.0

echo "PATH: $PATH"
which Rscript
Rscript --version
