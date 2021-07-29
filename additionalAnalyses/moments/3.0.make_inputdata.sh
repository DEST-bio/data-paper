#!/usr/bin/env bash
#

module load gcc/7.1.0 openmpi/3.1.4 R/3.6.3 anaconda/2020.11-py3.8

R_script=$1
Pops=$2
Out_dir=$3
Meta_dir=$4
DEST_gds_dir=$5

echo "began at"  `date`

cat $R_script | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d"

### format data for moments
Rscript \
$R_script \
${SLURM_ARRAY_TASK_ID} \
$Pops \
$Out_dir \
$Meta_dir \
$DEST_gds_dir