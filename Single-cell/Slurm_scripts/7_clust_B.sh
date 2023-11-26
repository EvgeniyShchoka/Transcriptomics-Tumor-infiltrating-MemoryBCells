#!/bin/sh -l

#SBATCH --job-name=7_clust_B        # Job name
#SBATCH --cpus-per-task=40         # Run on a single CPU
#SBATCH --mem=256G                 # Job memory request
#SBATCH --partition=short           # Time limit hrs:min:sec
#SBATCH --output=7_clust_B.%j.log   # Standard output and error log

module load R
Rscript --vanilla ~/B_Memory_master/Master_Diploma_Private/Output/Single-cell/R_scripts/7_clusterization_of_B_cells.R