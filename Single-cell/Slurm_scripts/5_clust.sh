#!/bin/sh -l

#SBATCH --job-name=5_clust        # Job name
#SBATCH --cpus-per-task=40         # Run on a single CPU
#SBATCH --mem=256G                 # Job memory request
#SBATCH --partition=medium           # Time limit hrs:min:sec
#SBATCH --output=5_clust.%j.log   # Standard output and error log

module load R
Rscript --vanilla ~/B_Memory_master/Master_Diploma_Private/Output/Single-cell/R_scripts/5_clusterization_of_all_cells.R