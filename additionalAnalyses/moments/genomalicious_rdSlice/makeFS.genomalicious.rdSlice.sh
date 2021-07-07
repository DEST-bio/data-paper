#!/usr/bin/env bash
#
#SBATCH -J makeFS # A single job name for the array
#SBATCH --ntasks-per-node=3 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-01:00  ### 10 minutes
#SBATCH --mem 24G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/makeFS.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/makeFS.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### run as: sbatch --array=1-$( cat /project/berglandlab/moments/pairs.csv | sed '1d' | wc -l ) /scratch/aob2x/data-paper/additionalAnalyses/moments/genomalicious_rdSlice/makeFS.genomalicious.rdSlice.sh
### run as: sbatch --array=1-4 /scratch/aob2x/data-paper/additionalAnalyses/moments/makeFS.sh

### sacct -j 23417985
### cat /scratch/aob2x/dest/slurmOutput/makeFS.23409129_1.out
module load gcc/7.1.0 openmpi/3.1.4 R/3.6.3

## SLURM_ARRAY_TASK_ID=2

Rscript /scratch/aob2x/data-paper/additionalAnalyses/moments/genomalicious_rdSlice/makeFS.genomalicious.rdSlice.R ${SLURM_ARRAY_TASK_ID}
#Rscript /scratch/aob2x/data-paper/additionalAnalyses/moments/makeFS.R 1
