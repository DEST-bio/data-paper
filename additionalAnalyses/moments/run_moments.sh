#!/usr/bin/env bash
#
#SBATCH -J makeFS # A single job name for the array
#SBATCH --ntasks-per-node=3 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-01:00  ### 10 minutes
#SBATCH --mem 24G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/run_moments.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/run_moments.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### run as: sbatch --array=1-$( cat /scratch/aob2x/data-paper/additionalAnalyses/moments/pairs.csv | sed '1d' | wc -l ) /scratch/aob2x/data-paper/additionalAnalyses/moments/run_moments.sh
### run as: sbatch --array=1-4 /scratch/aob2x/data-paper/additionalAnalyses/moments/makeFS.sh

### sacct -j 23372806
### cat /scratch/aob2x/dest/slurmOutput/makeFS.22880380_2.out
module load gcc/7.1.0 openmpi/3.1.4 R/3.6.3 anaconda/2020.11-py3.8

## SLURM_ARRAY_TASK_ID=3

echo "began at"  `date`

cat /scratch/aob2x/data-paper/additionalAnalyses/moments/pairs.csv | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d"

### format data for moments
  Rscript /scratch/aob2x/data-paper/additionalAnalyses/moments/makeSFS_data.R ${SLURM_ARRAY_TASK_ID}

### define parameters
  metadata=/scratch/aob2x/moments_general/input/${SLURM_ARRAY_TASK_ID}.meta
  cat ${metadata}

  Pair=$( cat $metadata  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $1 }' )
  SFS=$( cat $metadata  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $2 }' )
  L=$( cat $metadata  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $3 }' )
  pop1_id=$( cat $metadata  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $4 }' )
  pop2_id=$( cat $metadata  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $5 }' )

  projection1=$( cat $metadata  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $6 }' )
  projection2=$( cat $metadata  | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $7 }' )

  echo "Now Processing" $Pair
  echo "Now Loading" $SFS "=> where" $Pair "SFS is located"
  echo $Pair "has an L parameter of" $L "bp"


### run moments
  source activate moments_kern

  cd /scratch/aob2x/moments_general/output

  python /scratch/aob2x/data-paper/additionalAnalyses/moments/genomalicious/moments_genom_test2.py \
  ${SFS} \
  $L \
  100 \
  $Pair \
  $pop1_id \
  $pop2_id \
  $projection1 \
  $projection2

  conda deactivate

### Print the time
  echo "ended at"  `date`
