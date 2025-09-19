#!/bin/bash
#SBATCH --account=p32397
#SBATCH --partition=short
#SBATCH --time=00:05:00
#SBATCH --job-name=test_serial
#SBATCH --output=test_serial_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=500M

# Load R module
module purge all
module load R/4.4.0

# Run the serial script
Rscript test_serial.R
