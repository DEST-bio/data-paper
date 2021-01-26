#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=60G
#SBATCH --time=5:00:00
#SBATCH --partition=largemem
#SBATCH --account=jcbnunez
#SBATCH --array=1-265

module load gcc/7.1.0  
module load openmpi/3.1.4
module load gdal
module load proj
module load goolf R/4.0.0

Rscript LOOCV_Cluster.r ${SLURM_ARRAY_TASK_ID}

