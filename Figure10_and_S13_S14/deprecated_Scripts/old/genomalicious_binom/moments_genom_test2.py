# This script will run the program moments by Jouganous, J., Long, W., Ragsdale, A. P., & Gravel, S. (2017)
# Written by Keric Lamb, UVA 2021
# ksl2za@virginia.edu 
# import packages that'll be used
import moments
from moments import Numerics
from moments import Integration
from moments import Spectrum
from moments import Misc
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import  pandas as pd

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

projection1= int(projection1)
projection2= int(projection2)

#opening output file to give column names
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
            str("-2LL_model")+'\t'+
            str("AIC")+'\n')
PMmod.close()

#read in data as file
dd = Misc.make_data_dict(fs_file) #reads in genomalicious SNP file

pop_id=[pop_id1,pop_id2]
projection=[projection1,projection2]

fs_folded = Spectrum.from_data_dict(dd, pop_ids=pop_id, projections=projection, polarized=False) #takes data dict and folds
ns = fs_folded.sample_sizes #gets sample sizes of dataset

#this code was used to import Alan's custom fs files. Now commented out as we've switched to genomalicious (R package) compiled SNP files
# data = moments.Spectrum.from_file(fs_file)
# fs_folded = data.fold()
# ns = fs_folded.sample_sizes

# Highlighted out, as pts_l is only relevant for DaDi-- moments doesn't use grid-point extrapolations
# These are the grid point settings will use for extrapolation.
# They're somewhat arbitrary, but should be big enough that they don't change LL if increased
# pts_l = [100,140,180]

#model with symmetric migration from Moments bitbucket
def split_mig_moments(params, ns, pop_ids=None):
    if pop_ids is not None and len(pop_ids) != 2:
        raise ValueError("pop_ids must be a list of two population IDs")
    nu1, nu2, T, m = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T, m=np.array([[0, m], [m, 0]]))
    fs.pop_ids = pop_ids
    return fs

func_moments = split_mig_moments

#boundaries
upper_bound = [10, 10, 5, 50]
lower_bound = [1e-5, 1e-5, 1e-5, 0]

#constants
mu = 2.8e-9 #from Keightley et al. 2014
L = int(L_file) #imported from system argument #2
g = 0.07692308 #equals 13 gen/year. Calculated based on biological intuition.

# run X optimizations from uniform-sampled starting params
for i in range(int(iterations)): #iterations is imported from sys. argument #1
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
        str(theta)+'\t'+
        str(ll_model)+'\t'+
        str(aic)+'\n')
    PMmod.close()

print("Moments finished running")