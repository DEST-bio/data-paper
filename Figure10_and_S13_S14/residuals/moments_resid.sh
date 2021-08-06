#!/usr/bin/bash
#
#SBATCH -J resids # A single job name for the array
#SBATCH -N 1 # on one node
#SBATCH -n 1
#SBATCH --ntasks=1 # one core
#SBATCH -t 2:00:00
#SBATCH --mem=40G
#SBATCH -o ~/DEST_playground/out/moments_resids.out
#SBATCH -e ~/DEST_playground/out/moments_resids.error
#SBATCH -p standard
#SBATCH -a berglandlab

### sbatch --array=1-3 moments_resid.sh

echo "began at"  `date`

#Load conda module
module purge
module load anaconda/2020.11-py3.8

#Activate dadi kernel
source activate dadi_env

moments=~/DEST_playground/moments_resids_theta.py

#Load the metadata object into memory
metadata=/home/ksl2za/DEST_playground/resid_meta.txt #Address to the metadata. 

#Pair_names=> the names of the two pools being considered. Separated by a predictable delimiter like "|" (do not use "_")!!!
#fs=> the address of the file containing the 2D SFS
#iterations=> number of times for the script to loop

#===> Want to fix SLURM_ARRAY_TASK_ID? this is useful for debugging <====
#=# SLURM_ARRAY_TASK_ID=120

#Mining the metadata file
Pair=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $2 }' )
fsfile=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $1 }' )
param=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $3 }' )
model_type=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $4 }' )
pop1=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $5 }' )
pop2=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $6 }' )
proj1=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $7 }' )
proj2=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $8 }' )
nu1=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $9 }' )
nu2=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $10 }' )
m12=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $11 }' )
Ts=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $12 }' )
theta=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $13 }' )
nu1B=$$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $14 }' )
nu2B=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $15 }' )
nu1F=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $16 }' )
nu2F=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $17 }' )

echo "Now Processing" $Pair
echo $Pair "is running" $iter "iterations"

#Run moments
#fs_file = sys.argv[1]
#Pair_name = sys.argv[2]
#param = sys.argv[3]
#model_type= sys.argv[4]
#pop_name1 = sys.argv[5]
#pop_name2 = sys.argv[6]
#pool_n1 = sys.argv[7]
#pool_n2 = sys.argv[8]
#nu1 = float(sys.argv[9])
#nu2 = float(sys.argv[10])
#m12 = float(sys.argv[11])
#Ts = float(sys.argv[12])
#theta = float(sys.argv[13])
#nu1B = float(sys.argv[14])
#nu2B = float(sys.argv[15])
#nu1F = float(sys.argv[16])
#nu2F = float(sys.argv[17])

cd ~/DEST_playground/
MPLBACKEND=Agg python $moments \
$fsfile \
$Pair \
$pop1 \
$pop2 \
$proj1 \
$proj2 \
$nu1 \
$nu2 \
$m12 \
$Ts \
$theta \
$nu1B \
$nu2B \
$nu1F \
$nu2F

#De-Activate moments kernel
conda deactivate

#Print the time
echo "ended at"  `date`