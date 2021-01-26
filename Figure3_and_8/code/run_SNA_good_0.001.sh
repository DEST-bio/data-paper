#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=60G
#SBATCH --time=8:00:00
#SBATCH --partition=standard
#SBATCH --account=jcbnunez

module load gcc/7.1.0  
module load openmpi/3.1.4
module load gdal
module load proj
module load htslib/1.9
module load bcftools/1.9
module load intel/18.0  
module load intelmpi/18.0
module load R/4.0.0

Rscript 2.1.AlleleFreq.SNAPE.goodsamps.0.001.R