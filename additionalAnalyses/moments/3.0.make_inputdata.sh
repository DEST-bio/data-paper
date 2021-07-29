#!/usr/bin/env bash
#

module load intel/18.0 intelmpi/18.0  
module load goolf/7.1.0_3.1.4  
module load gdal proj R/4.0.0
module load anaconda/2020.11-py3.8


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