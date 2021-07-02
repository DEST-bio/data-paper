
# Set up environemnt
#load dependencies
module load anaconda/2020.11-py3.8

#Step 1. Create conda environemt
#This needs to be done  only once

conda create \
-n moments_kern \
python=3.8 \
dadi \
ipykernel \
-c conda-forge

# Activates conda environemnt
source activate moments_kern

#Installs moments
conda install moments -c bioconda
