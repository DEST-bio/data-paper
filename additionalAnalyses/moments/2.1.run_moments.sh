#!/usr/bin/env bash
#
#SBATCH -J makeFS # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-06:00  ### 10 minutes
#SBATCH --mem 18G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/makeFS.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/makeFS.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch --array=1-$( wc -l /project/berglandlab/moments/moments_jobs.delim ) /scratch/aob2x/data-paper/additionalAnalyses/moments/2.1.run_moments.sh
### sbatch --array=1-126 /scratch/aob2x/data-paper/additionalAnalyses/moments/2.1.run_moments.sh
### sacct -j 22847954

echo "began at"  `date`

#Load conda module
module load anaconda/2020.11-py3.8

#Activate moments kernell
source activate moments_kern

#Load the metadata object into memory
metadata=/project/berglandlab/moments/moments_jobs.delim #Address to the metadata. What is this? see below:

#The metadata file is as follows:
# A file with 3 columns and a header
# Separated by tabs
#Pair_names \t address_to_fs \t L
#At1|At2	   /file/At12.fs   1200000

#Pair_names=> the names of the two pools being considered. Separated by a predictable delimiter like "|" (do not use "_")!!!
#address_to_fs=> the address of the file containing the 2D SFS
#L => The number of sites from which the SFS was derived. IMPORTANT: This number must be an integer. No scientific notation allowed.

#===> Want to fix SLURM_ARRAY_TASK_ID? this is useful for debugging <====
#=# SLURM_ARRAY_TASK_ID=120

#Mining the metadata file

Pair=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $1 }' )
SFS=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $2 }' )
L=$( cat $metadata  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk -F "\t" '{ print $3 }' )

echo "Now Processing" $Pair
echo "Now Loading" $SFS "=> where" $Pair "SFS is located"
echo $Pair "has an L parameter of" $L "bp"

#Run Moments
# this script takes 4 arguments
#fs_file = sys.argv[1] ==> location of SFS i.e. -> $SFS
#L_file = sys.argv[2] ==> L param i.e. -> $L
#iterations = sys.argv[3] ==> number of runs ... a number
#pair_name = sys.argv[4] ==> name of the pair -> $Pair

cd /project/berglandlab/moments/moments_output
python /scratch/aob2x/data-paper/additionalAnalyses/moments/2.2.MomentsCode.py \
$SFS \
$L \
100 \
$Pair


#De-Activate moments kernell
conda deactivate

#Print the time
echo "ended at"  `date`
