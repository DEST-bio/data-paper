# Instructions for demographic inference using moments and the DEST pool-seq data

# Description
This analysis will use the demographic estimator moments to analyze the demographic history and population genetic parameters from populations of European drosophila part of the DEST dateset. In particular, we focus on the E/W split described in the main text. We have implemented a simple "split-w/-migration" model to estimate divergence time.

# Scripts in this analysis
Here is a summary of the scripts presented in this analysis. **Important**: all analyses were done using SLURM headers which are customized to our cluster. You will have to modify these headers so that they may run in your own cluster.

0. **set_up_moments.RunManually.sh**: this script, meant to be run manually, install the program moments in your machine. You may choose to install it your own way.
1. **generatePairs.R**: This script draws from the metadata of the DEST repository (https://github.com/DEST-bio/DEST_freeze1) to generate all possible comparisons between populations in the European clusters.
2. **getL.R**: This script also draws from the DEST repo and its used to estimate the number of bases sequenced used to estimate the SFS.
3. **make_inputdata.sh**: This script implements an array slurm job to generate the input files for moments. This script has an R companion (3.0.make_inputdata.sh + 3.1.makeSFS_data.R)
4. **run_moments.sh**: This is the wrapper script which runs the program moments. This wrapper calls on two possible versions of the python moments script: 4.1.moments_param_theta.py + 4.2.moments_bounds.py 
5. **collect_all_moments.R**
6. **plot_all_moments.R**

# Script 0. Install moments
Here we provide you with an example to install the python libraries needed to run moments -- Install moments (`data-paper/additionalAnalyses/moments/0.set_up_moments.RunManually.sh`)

# Script 1. Generate Pairs
Because the East and West cluster are composed by  populations geographically distributed across west and east Europe, we want to estimate diverge-w/-migration accounting for this spatial separation. For example, we hope to estimate and compare parameters from populations coming from East and West Austria (the boundary region) but also, Spain vs Ukraine (two spatially distant populations within their clades). This script will generate 1000 random pairs among estern and western flies for the analysis. -- Generate ~1000 pairs of populations between E and W European clusters (`data-paper/additionalAnalyses/moments/1.generatePairs.R`).

### Note: Other repos are required
This script draws from the metadata file from the dest dataset "DEST_freeze1/populationInfo/samps_10Nov2020.csv" to create the population pairs. As such, you must download the git repo of the DEST datset.

```bash
cd <../"where I want my repo to be">
git clone https://github.com/DEST-bio/DEST_freeze1.git
```
Make sure you modify the script such that the file is being pulled from the right location. 
```bash
  samps <- fread("<repo>/DEST_freeze1/populationInfo/samps_10Nov2020.csv")
```

### What is the script output?
This script outputs various files:
```R
"./pairs_between.txt"
"./pairs_within.txt"
"./pairs_all.txt"
``` 
These files contained different types of comparisons within and between populations in Europe. However, for the purposes of replicating our main text analysis. The file you want is **"./pairs_all.txt"**. This file will instruct the code to do all possible--and relevant--comparisons, 

# Script 2. Get L for moments calculation
In order to estimate divergence time and other parameters, it is needed to have an estimate of how many bases were lest unmasked after filtration. -- Get median L for moments calculations (`data-paper/additionalAnalyses/moments/2.getL.R`). Same as above, make sure you modify the script such that the file is being pulled from the right location. 

# Script 3. Make input data for moments
Because we intend to run 1000+ comparisons, models, and demographic estimation, we will run the program moments using an array configuration in a HPC environment. This require us to use an array job. Example:

```bash
sbatch --array=1-$( cat /<workfolder>/pairs_all.all.csv | sed '1d' | wc -l ) \
/<workfolder>/3.0.make_inputdata.sh
```
In this case, we are running an array with as many jobs as pairs in the generate pairs output.

### What is the script doing? 

This script implements the R package *genomalicious* (https://github.com/j-a-thia/genomalicious) to discretize allele frequencies from Pool-Seq data. This is a fundamental requirement for SFS estimation and running moments. Similar to our other scripts above, this script depends on data from the DEST repo, as well as data from the DEST pipeline output. So make sure to update the script and link it to the proper metadata file:

```bash
#data from the DEST repo:
samps <- fread("<repo>/DEST_freeze1/populationInfo/samps_10Nov2020.csv")
#data from the DEST pipeline output. You will have to run all the DEST pipeline to reproduce these results!
genofile <- seqOpen(paste("/project/berglandlab/DEST/gds/dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.gds", sep=""))
#or
genofile <- seqOpen(paste("/project/berglandlab/DEST/gds/dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.gds", sep=""))
```
You will also need to install R packages

```R
library(data.table)
library(SeqArray)
library(foreach)
library(sp)
library(doMC)
library(genomalicious)
```

### What are the script options?
Aside from the *SLURM_ARRAY_TASK_ID* input, which is derived from the array job itself, we originally designed this script to run comparisons for: all samples, between-, or within-clusters. These options are, respectively, *all_pops, between,* and *within*. As mentioned above, the only option which is probably relevant to you is ***all_pops***. Example:
```bash
Rscript /data-paper/additionalAnalyses/moments/3.1.makeSFS_data.R \
${SLURM_ARRAY_TASK_ID} \
all_pops
```

### What is the script output?
The script will create a directory, and will output all input files for moments there. 

# Script 4. Run *moments*
In this step, we will run the script which calls on *moments* to do demographic inference. This code also runs as an array:

```bash
sbatch --array=1-$( cat /scratch/aob2x/data-paper/additionalAnalyses/moments/pairs_between.csv | sed '1d' | wc -l ) \
./4.0.run_moments.sh \
pairs_all.txt \
4.1.moments_bounds.py \
50
/project/berglandlab/moments/moments_input
```

### Whats inside the *moments* wrapper
The script 5.run_moments_between.sh runs a wrapper of the python code for moments. However, because we are running the program as an array the script provides python with different data for each population pair. For example, we first load all the packages needed to activate python and run the array
```
#Load needed modules
module load gcc/7.1.0 
openmpi/3.1.4 R/3.6.3 
anaconda/2020.11-py3.8
```

### Provide the metadata, analysis type, and iterations num
For this script to work it is important to provide two user defined inputs, the metadata file, and the type of analysis to run.
```bash
sbatch --array=1-$( cat /scratch/aob2x/data-paper/additionalAnalyses/moments/pairs_between.csv | sed '1d' | wc -l ) \
./4.0.run_moments.sh \
pairs_all.txt \ <=== Metadata file
4.1.moments_bounds.py \ <=== Analysis Type
50 \ <=== Fifty iterations each run
/project/berglandlab/moments/moments_input <=== Input data
```
* The metadata can be: "./pairs_all.txt", "./pairs_between.txt", "./pairs_within.txt"
* The analyses are the actual python scripts and can be: **"4.1.moments_bounds.py"** or **"4.2.moments_param_theta.py"**. The first, *moments_bounds.py* is close to the *default* way to run moments. And, *moments_param_theta.py* is our modified script which parametrized theta (see main manuscript).

### Importing data into the moments script
This wrapper will import data from two sources. The metadata file, and the output of 3.0.make_inputdata.sh.

```bash
## User provided variables
metadata=$1
moments=$2
iterations=$3
input_folder=$4

#generate input variables for moments
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
```

## 5.3. Run python--moments 
Finally, the wrapper runs the python script for moments. This is run using a custom kernel. The custom kernel was generated in step zero (see above). The variable "$moments" is the version of analysis that 

```bash
### run moments
  source activate moments_kern
  
  mkdir output
  
  cd ./output

  python $moments \
  ${SFS} \
  $L \
  $iterations \
  $Pair \
  $pop1_id \
  $pop2_id \
  $pool_n1 \
  $pool_n2

  conda deactivate

```
Notice that:
```bash
  ${SFS} \ ==> this is the SFS
  $L \ ==> this is the number of sequenced sites used
  $iterations \ ==> This is the number of iterations moments runs for
  $Pair \ ==> This is the name of the pair
  $pop1_id \ ==> This is the id of pop1
  $pop2_id \ ==> This is the id of pop2 
  $pool_n1 \ ==> This is the number of inds. pooled in pop1
  $pool_n2 ==> This is the number of inds. pooled in pop2
```
### Moments --> Dependencies

The script uses standard libraries to run:
```python
import moments
from moments import Numerics
from moments import Integration
from moments import Spectrum
from moments import Misc
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
```
### Moments --> Import parameters
Inside the python script these are passed as sys.arg variables
```
#define sys args
#DEBUGGED AS FOR LOOPS
fs_file = sys.argv[1]
L_file = sys.argv[2]
iterations = sys.argv[3]
Pair_name = sys.argv[4]
pop_id1 = sys.argv[5]
pop_id2 = sys.argv[6]
projection1 = sys.argv[7]
projection2 = sys.argv[8]
```
### Moments --> The output file
The output of moments will be capture in a table with the following data. 
* Pair_name: The names of the population pair
* fs_name: prints de delim file which was used in the model
* L: print the L parameter
* pop1_size: This the effective Ne of pop1 (model output)
* pop2_size: This the effective Ne of pop2 (model output)
* divergence_time: This is the time in years that pop1 and pop2 diverged (model output)
* mig_pop1: This the number of individuals in pop1 that derived from pop2  (model output)
* mig_pop2: This the number of individuals in pop2 that derived from pop1  (model output)
* theta: Scales parameter used to scale the SFS (4NeμL)
* nu1: Raw model output that will scale the ancestral pop. size to the size of the new pop1
* nu2: Raw model output that will scale the ancestral pop. size to the size of the new pop2
* Ts: raw output for divergence time
* m12: raw output for migration rate
* -2LL_model: Composite log likelihood. closer to zero is better 
* AIC: information criterion.  Weights the model assuming 4 parameters. closer to zero is better.
```python
PMmod=open('%s_output.txt' % Pair_name,'a')
PMmod.write(
            str("Pair_name")+'\t'+ #print pair name
            str("fs_name")+'\t'+ #double checking fs_lines[y] is working as I want it to
            str("L")+'\t'+ #double checking L is working as I want it to
            str("pop1_size")+'\t'+ #nu1
            str("pop2_size")+'\t'+ #nu2
            str("divergence_time")+'\t'+ #divergence T
            str("mig_pop1")+'\t'+ #Migration ij
            str("mig_pop2")+'\t'+ #Migration ji
            str("theta")+'\t'+
            str("nu1")+'\t'+
            str("nu2")+'\t'+
            str("Ts")+'\t'+
            str("m12")+'\t'+
            str("-2LL_model")+'\t'+
            str("AIC")+'\n')
PMmod.close()

```
### Moments --> Import data
Make a data dictionary from the SFS file 
```
dd = Misc.make_data_dict(fs_file) #reads in genomalicious SNP file
```
Load in information about the data. Here, the variable projection1 means the pool sizes.
```
pop_id=[pop_id1,pop_id2]
projection=[projection1,projection2]
```
The model is loading the SFS and folding it. This is done to simplify the calculation in absence of ancestral states.
```
fs_folded = Spectrum.from_data_dict(dd, pop_ids=pop_id, projections=projection, polarized=False) #takes data dict and folds
ns = fs_folded.sample_sizes #gets sample sizes of dataset
```
### Moments --> The guts of the mode "split-w/-migration"
We are defining a function which defines the model itself
```python
def split_mig_moments(params, ns, pop_ids=None):
    if pop_ids is not None and len(pop_ids) != 2:
        raise ValueError("pop_ids must be a list of two population IDs")
    nu1, nu2, T, m = params

#"moments.LinearSystem_1D.steady_state_1D" is a premade function from the moments bitbucket --> The basic behaviour of the model is to model the pop1 and pop2 as a single evolving populations under a steady state 
#moments.Spectrum --> creates a spectrum object
#moments.Manips.split_1D_to_2D --> is splitting the population into 2
#integrate --> this step defines the conditions of the split. 
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T, m=np.array([[0, m], [m, 0]]))
    fs.pop_ids = pop_ids
    return fs

func_moments = split_mig_moments
```
The model is a split with migration type. Which means that we are assuming that one population splint into two populations, some time in the past, but migration among the two subdivisions continues at a symmetric rate. If the user attempt to fit the model for more than two pops an error will generate.
The model has 4 parameters (nu1, nu2, T, m) which will be estimated by the model. The nu1, nu2 params are the fractions of the ancestral population (Nref) after the split. These fractions do not need to add up to 1. 

### Moments --> Setting boundary behavior
We will set the boundaries from which the model can draw parameters and also constrain the model optimization. Limitations are in place to prevent pathological runaway parameter behavior which tends to converge toward zero and infinity.  This is an example from: "4.1.moments_bounds.py". 
```python
#boundaries for nu1, nu2, Ts, m12
upper_bound = [10, 10, 5, 50]
lower_bound = [1e-2, 1e-2, 1e-2, 0]
```
In the case of 4.2.moments_param_theta.py there is an aditional parameter for theta:

```python
#boundaries for nu1, nu2, Ts, m12, theta1
upper_bound = [10, 10, 5, 50, 5e5]
lower_bound = [1e-5, 1e-5, 1e-5, 0, 3.75e5]
```

### Moments --> Providing conversion constants
We provide biological constants for μ (mutation rate), L, and g (generation time). These are obtained from the literature and biological intuition.
```
#constants
mu = 2.8e-9 #from Keightley et al. 2014
L = int(L_file) #imported from system argument #2
g = 0.07692308 #equals 13 gen/year. Calculated based on biological intuition.
```
### Moments --> Running the model
Finally, the script runs the model inside a loop. The loop will run the number of times as given in argument #3 in the wrapper (50 in this case).
```python
# run X optimizations from uniform-sampled starting params
for i in range(int(iterations)): #iterations is imported from sys. argument named "iterations"
    print("starting optimization "+str(i))

    #Start the run by picking random parameters from a uniform distribution.
    #Parameters are set above " nu1, nu2, T, m"
    #The number "4" in range(4) comes from the number of parameters. Change in needed.
    popt=[np.random.uniform(lower_bound[x],upper_bound[x]) for x in range(4)]

    #This is the optimization step for moments.
    #popt is the prior.
    #fs folded is a tranform SFS by folding it. The original SFS loaded in sys.arg #1 is quasi-folded. i.e. polarized to reference genome.
    #Folding it is a must, because reference is not 100% ancestral
    popt=moments.Inference.optimize_log(popt, fs_folded, func_moments,
                                        lower_bound=lower_bound, upper_bound=upper_bound,
                                        verbose=False, maxiter=100,
                                        )

    #This number is 4. i.e., count parameters. there is an opportunity to streamline the code by propagating this from the beginning.
    params = len(["nu1", "nu2", "Ts", "m"]) #for use in AIC calculation

    #This is the moments function.
    model = func_moments(popt, ns)

    #Calculate log likelihood of the model fit
    ll_model=moments.Inference.ll_multinom(model, fs_folded)
    #Now calculate AIC of model fit
    aic = 2*params - 2*ll_model
    print('Maximum log composite likelihood: {0}'.format(ll_model))

    #Now estimate theta from model fit
    theta = moments.Inference.optimal_sfs_scaling(model, fs_folded)
    #Now calculate Ts from Model fit
    divergence_time = 2*(theta/(4*mu*L))*popt[2]*g #calculates divergence time in years

    #Now calculate Migration rate (fraction of migrants that move between pops)
    Mij = popt[3]/(2*(theta/(4*mu*L))) #actual migration rate
    #Below is an old code which had a typo. Keric has since updated it. kept for record keeping. #Jcbn Jun18,2021
    #mig_pop1 = Mij*(2*popt[0]) #number of individuals going i to j
    #mig_pop2 = Mij*(2*popt[1]) #number of individuals going j to i

    #Now we are estimated the nominal migration rate based on Mij
    mig_pop1 = Mij*(popt[0]*(theta/(4*mu*L))) #number of individuals going i to j: migrants = Mij*nu1*Nref
    mig_pop2 = Mij*(popt[1]*(theta/(4*mu*L))) #number of individuals going j to i: migrants = Mij*nu2*Nref

    #Now estimate population size
    pop1_size = popt[0]*(theta/(4*mu*L)) #pop1 size
    pop2_size = popt[1]*(theta/(4*mu*L)) #pop2 size

    #Open the output file
    PMmod=open('%s_output.txt' % Pair_name,'a')

    #Dumping output ot outfile
    PMmod.write(
        str(Pair_name)+'\t'+ #print pair name
        str(fs_file)+'\t'+ #double checking fs is the right one
        str(L)+'\t'+ #double checking L is working as desired
        str(pop1_size)+'\t'+ #nu1
        str(pop2_size)+'\t'+ #nu2
        str(divergence_time)+'\t'+ #divergence T
        str(mig_pop1)+'\t'+ #Migration ij
        str(mig_pop2)+'\t'+ #Migration ji
        str(theta_model)+'\t'+ #theta as calculated by the model
        str(theta_param)+'\t'+ #theta as calculated as a parameter. should equal theta model, just doing sanity check
        str(nu1)+'\t'+ #raw parameter output
        str(nu2)+'\t'+ #raw parameter output
        str(Ts)+'\t'+ #raw parameter output
        str(m12)+'\t'+ #raw parameter output
        str(ll_model)+'\t'+
        str(aic)+'\n')
    PMmod.close()

print("Moments finished running")
```
# Script 6. Collate all outputs 
Since we are running the moments script over multiple pairs we need to centralize all the  outputs into a single file for data mining. This is done here --  Collect moments data (`data-paper/additionalAnalyses/moments/6.collect_all_moments.R`)


# Script 7. Plot and misc.
We provide scripts for plotting and misc. analysis as shown in the paper. Although we recommend the end user to do this themselves for other analyses. -- Plot moments data (`data-paper/additionalAnalyses/moments/7.plot_all_moments.R`)
