#!/bin/sh -l

#SBATCH --job-name=8_batch_B        # Job name
#SBATCH --cpus-per-task=40         # Run on a single CPU
#SBATCH --mem=256G                 # Job memory request
#SBATCH --partition=short           # Time limit hrs:min:sec
#SBATCH --output=8_batch_B.%j.log   # Standard output and error log

module load R
Rscript --vanilla ~/B_Memory_master/Master_Diploma_Private/Output/Single-cell/R_scripts/8_batch_correction_of_B_cells_and_clusterization.R