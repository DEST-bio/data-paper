#!/usr/bin/env bash
#

### run as: sbatch --array=1-$( cat /scratch/aob2x/data-paper/additionalAnalyses/moments/pairs_all.csv | sed '1d' | wc -l ) /scratch/aob2x/data-paper/additionalAnalyses/moments/run_moments.sh

### between
### sacct -j 23454077
### cat /scratch/aob2x/dest/slurmOutput/run_moments.23453934_5.err


module load gcc/7.1.0 openmpi/3.1.4 R/3.6.3 anaconda/2020.11-py3.8

## User defined variables

metadata=$1
moments=$2
iterations=$3
input_folder=$4

echo "began at"  `date`

### pull parameters from the metadata flie

  pop1_id=$( cat $metadata | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d"  | awk -F " " '{ print $5 }' )
  echo $pop1_id
  
  pop2_id=$( cat $metadata | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d"  | awk -F " " '{ print $6 }' )
  echo $pop2_id

  Pair=$(echo ${pop1_id}.${pop2_id} )
  echo $Pair
	
  SFS_method=$( cat $metadata | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d"  | awk -F " " '{ print $3 }' )
  echo $SFS_method

  Caller=$( cat $metadata | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d"  | awk -F " " '{ print $2 }' )
  echo $Caller

  Demo=$( cat $metadata | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d"  | awk -F " " '{ print $11 }' )
  echo $Demo

### Sample to be analyzed
head $input_folder/$Caller.$SFS_method.$Pair.$Demo.meta
head $input_folder/$Caller.$SFS_method.$Pair.$Demo.delim

### Prepare additional metadata
  SFS=$input_folder/$Caller.$SFS_method.$Pair.$Demo.delim

  L=$( cat $input_folder/$Caller.$SFS_method.$Pair.$Demo.meta   | awk -F "\t" '{ print $3 }' )

  pool_n1=$( cat $input_folder/$Caller.$SFS_method.$Pair.$Demo.meta  | awk -F "\t" '{ print $6 }' )
  pool_n2=$( cat $input_folder/$Caller.$SFS_method.$Pair.$Demo.meta  | awk -F "\t" '{ print $7 }' )

### Finally a sanity check
  echo "Now Processing" $Pair
  echo "Now Loading" $SFS "=> where" $Pair "SFS is located"
  echo $Pair "has an L parameter of" $L "bp"


### run moments
  source activate moments_kern
  
  mkdir output
  
  cd ./output

  python $moments \
  ${SFS} \
  $L \
  $iterations \
  $Caller.$SFS_method.$Pair.$Demo \
  $pop1_id \
  $pop2_id \
  $pool_n1 \
  $pool_n2

  conda deactivate

### Print the time
  echo "ended at"  `date`
