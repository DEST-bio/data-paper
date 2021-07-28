#!/usr/bin/env bash
#
#SBATCH -J makeFS # A single job name for the array
#SBATCH --ntasks-per-node=3 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-01:00  ### 8 hours / job
#SBATCH --mem 27G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/run_moments.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/run_moments.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### run as: sbatch --array=1-$( cat /scratch/aob2x/data-paper/additionalAnalyses/moments/pairs_all.all.csv | sed '1d' | wc -l ) /scratch/aob2x/data-paper/additionalAnalyses/moments/make_inputdata.sh

### within
### sacct -j 23745161
### cat /scratch/aob2x/dest/slurmOutput/run_moments.23453934_5.err

module load gcc/7.1.0 openmpi/3.1.4 R/3.6.3 anaconda/2020.11-py3.8

## SLURM_ARRAY_TASK_ID=5

echo "began at"  `date`

cat /scratch/aob2x/data-paper/additionalAnalyses/moments/pairs_all.txt | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d"

### format data for moments
  Rscript /scratch/aob2x/data-paper/additionalAnalyses/moments/3.1.makeSFS_data.R ${SLURM_ARRAY_TASK_ID} all_pops
