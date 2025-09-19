#!/bin/bash
#SBATCH --account=p32397
#SBATCH --partition=short
#SBATCH --time=00:10:00
#SBATCH --job-name=test_multicore
#SBATCH --output=test_multicore_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=500M

module purge all
module load R/4.4.0

echo "Running minimal multicore future test..."
Rscript test_future_multicore.R
