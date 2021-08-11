#!/usr/bin/bash
#
#SBATCH -J resids # A single job name for the array
#SBATCH -N 1 # on one node
#SBATCH -n 1
#SBATCH --ntasks=1 # one core
#SBATCH -t 0:15:00
#SBATCH --mem=8G
#SBATCH -o moments_resids.out
#SBATCH -e moments_resids.error
#SBATCH -p standard
#SBATCH -a berglandlab

### sbatch --array=1-$( sed '1d' /project/berglandlab/moments/resid_meta.delim | wc -l ) /scratch/aob2x/data-paper/Figure10_and_S13_S14/residuals/moments_resid.sh

echo "began at"  `date`

#Load conda module
module purge
module load anaconda/2020.11-py3.8

#Activate dadi kernel
source activate moments_kern

moments=/scratch/aob2x/data-paper/Figure10_and_S13_S14/residuals/resids_all.v2.py

#Load the metadata object into memory
metadata=/project/berglandlab/moments/resid_meta.delim #Address to the metadata.

#Pair_names=> the names of the two pools being considered. Separated by a predictable delimiter like "|" (do not use "_")!!!
#fs=> the address of the file containing the 2D SFS
#iterations=> number of times for the script to loop

#===> Want to fix SLURM_ARRAY_TASK_ID? this is useful for debugging <====
#=# SLURM_ARRAY_TASK_ID=1

#Mining the metadata file
fsfile=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $1 }' )
Pair=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $2 }' )
param=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $3 }' )
model_type=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $4 }' )
model_sym=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $5 }' )
pop1=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $6 }' )
pop2=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $7 }' )
proj1=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $8 }' )
proj2=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $9 }' )
nu1=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $10 }' )
nu2=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $11 }' )
m12=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $12 }' )
m21=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $13 }' )
Ts=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $14 }' )
theta=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $15 }' )
nu1B=$$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $16 }' )
nu2B=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $17 }' )
nu1F=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $18 }' )
nu2F=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $19 }' )

echo "Now Processing" $Pair
echo $Pair "is running" $iter "iterations"

#Run moments
# fs_file = sys.argv[1]
# Pair_name = sys.argv[2]
# param = sys.argv[3]
# model_type= sys.argv[4]
# model_sym= sys.argv[5]
# pop_name1 = sys.argv[6]
# pop_name2 = sys.argv[7]
# pool_n1 = sys.argv[8]
# pool_n2 = sys.argv[9]
# nu1 = float(sys.argv[10])
# nu2 = float(sys.argv[11])
# m12 = float(sys.argv[12])
# m21 = float(sys.argv[13])
# Ts = float(sys.argv[14])
# theta = float(sys.argv[15])
# nu1B = float(sys.argv[16])
# nu2B = float(sys.argv[17])
# nu1F = float(sys.argv[18])
# nu2F = float(sys.argv[19])

#edit cd to desired location
#MPLBACKEND=Agg deals with matplotlib calling a display, which won't work in slurm
cd /project/berglandlab/moments/resids

MPLBACKEND=Agg python $moments \
$fsfile \
$Pair \
$param \
$model_type \
$model_sym \
$pop1 \
$pop2 \
$proj1 \
$proj2 \
$nu1 \
$nu2 \
$m12 \
$m21 \
$Ts \
$theta \
$nu1B \
$nu2B \
$nu1F \
$nu2F &

#De-Activate moments kernel
conda deactivate

#Print the time
echo "ended at"  `date`
