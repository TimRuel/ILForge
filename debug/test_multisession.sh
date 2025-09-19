#!/bin/bash
#SBATCH --account=p32397
#SBATCH --partition=short
#SBATCH --job-name=parallel_test
#SBATCH --mail-type=ALL
#SBATCH --mail-user=timothyruel2024@u.northwestern.edu
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err
#SBATCH --array=1-10             # 10 independent replications
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64        # Number of cores per task
#SBATCH --nodes=1
#SBATCH --ntasks=1 
#SBATCH --time=01:00:00
#SBATCH --mem=4G

# Load R module (adjust if needed)
module load R/4.4.0

# Set environment variables for multisession
export FUTURE_AVAILABLE_CORES=$SLURM_CPUS_PER_TASK

# Run R script with the array index as argument
Rscript test_multisession.R $SLURM_ARRAY_TASK_ID
