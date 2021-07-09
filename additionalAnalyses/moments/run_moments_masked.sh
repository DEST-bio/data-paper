#!/usr/bin/env bash
#
#SBATCH -J makeFS # A single job name for the array
#SBATCH --ntasks-per-node=3 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-08:00  ### 8 hours / job
#SBATCH --mem 27G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/run_moments.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/run_moments.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### cat /scratch/aob2x/data-paper/additionalAnalyses/moments/pairs_all.csv | grep -v "median_1sd" > /scratch/aob2x/data-paper/additionalAnalyses/moments/pairs_all.all.csv

### run as: sbatch --array=1-$( cat /scratch/aob2x/data-paper/additionalAnalyses/moments/pairs_all_all.csv | sed '1d' | wc -l ) /scratch/aob2x/data-paper/additionalAnalyses/moments/run_moments_masked.sh
### run as: sbatch --array=1-4 /scratch/aob2x/data-paper/additionalAnalyses/moments/run_moments_masked.sh

### between
### sacct -j 23528648
### cat /scratch/aob2x/dest/slurmOutput/run_moments.23528648_1.err


module load gcc/7.1.0 openmpi/3.1.4 R/3.6.3 anaconda/2020.11-py3.8

## SLURM_ARRAY_TASK_ID=5

echo "began at"  `date`

cat /scratch/aob2x/data-paper/additionalAnalyses/moments/pairs_all.all.csv | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d"

### format data for moments
  Rscript /scratch/aob2x/data-paper/additionalAnalyses/moments/makeSFS_data.R ${SLURM_ARRAY_TASK_ID} all_all

### define parameters
  metadata=/scratch/aob2x/moments_general/input_masked/${SLURM_ARRAY_TASK_ID}.meta
  cat ${metadata}

  Pair=$( cat $metadata  | awk -F "\t" '{ print $1 }' )
  SFS=$( cat $metadata  | awk -F "\t" '{ print $2 }' )
  L=$( cat $metadata  | awk -F "\t" '{ print $3 }' )
  pop1_id=$( cat $metadata  | awk -F "\t" '{ print $4 }' )
  pop2_id=$( cat $metadata  | awk -F "\t" '{ print $5 }' )

  projection1=$( cat $metadata  | awk -F "\t" '{ print $6 }' )
  projection2=$( cat $metadata  | awk -F "\t" '{ print $7 }' )

  echo "Now Processing" $Pair
  echo "Now Loading" $SFS "=> where" $Pair "SFS is located"
  echo $Pair "has an L parameter of" $L "bp"


### run moments
  source activate moments_kern

  cd /scratch/aob2x/moments_general/output_masked

  python /scratch/aob2x/data-paper/additionalAnalyses/moments/moments_genom_test2_mask.py \
  ${SFS} \
  $L \
  50 \
  $Pair \
  $pop1_id \
  $pop2_id \
  $projection1 \
  $projection2

  python /scratch/aob2x/data-paper/additionalAnalyses/moments/moments_genom_singlePop_mask.py \
  ${SFS} \
  $L \
  ${Pair} \
  $pop1_id \
  $pop2_id \
  $projection1 \
  $projection2

  conda deactivate

### Print the time
  rm $SFS
  rm $metadata
  echo "ended at"  `date`
