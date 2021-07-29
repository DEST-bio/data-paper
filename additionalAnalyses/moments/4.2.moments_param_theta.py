# This script will run the program moments by Jouganous, J., Long, W., Ragsdale, A. P., & Gravel, S. (2017)
# This script runs theta as a parameter of the model
# Written by Keric Lamb, UVA 2021
# ksl2za@virginia.edu
# import packages that'll be used
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import moments
from moments import Numerics
from moments import Integration
from moments import Spectrum
import dadi
from dadi import Misc

#define sys args
#DEBUGGED AS FOR LOOPS
fs_file = sys.argv[1]
L_file = sys.argv[2]
iterations = sys.argv[3]
Pair_name = sys.argv[4]
pop_id1 = sys.argv[5]
pop_id2 = sys.argv[6]
projection1 = sys.argv[8]
projection2 = sys.argv[7]

projection1= int(projection1)
projection2= int(projection2)

#opening output file to give column names
PMmod=open('./outputs_param/%s_output.txt' % Pair_name,'w')
PMmod.write(
            str("Pair_name")+'\t'+ #print pair name
            str("fs_name")+'\t'+ #double checking fs_lines[y] is working as I want it to
            str("L")+'\t'+ #double checking L is working as I want it to
            str("pop1_size")+'\t'+ #nu1
            str("pop2_size")+'\t'+ #nu2
            str("divergence_time")+'\t'+ #divergence T
            str("mig_pop1")+'\t'+ #Migration ij
            str("mig_pop2")+'\t'+ #Migration ji
            str("theta_model")+'\t'+ #sanity check, should equal ~1
            str("theta_param")+'\t'+ 
            str("nu1")+'\t'+
            str("nu2")+'\t'+
            str("Ts")+'\t'+
            str("m12")+'\t'+
            str("fs_sanitycheck")+'\t'+ #sanity checking that 
            str("-2LL_model")+'\t'+
            str("AIC")+'\n')
PMmod.close()

#read in data as file
dd = dadi.Misc.make_data_dict(fs_file) #reads in genomalicious SNP file

pop_id=[pop_id1,pop_id2]
projection=[projection1,projection2]

fs_folded = Spectrum.from_data_dict(dd, pop_ids=pop_id, projections=projection, polarized=False) #takes data dict and folds
ns = fs_folded.sample_sizes #gets sample sizes of dataset
S = fs_folded.S()
#fs_folded.mask[:1,:] = True
#fs_folded.mask[ :,:1] = True

#this code was used to import Alan's custom fs files. Now commented out as we've switched to genomalicious (R package) compiled SNP files
# data = moments.Spectrum.from_file(fs_file)
# fs_folded = data.fold()
# ns = fs_folded.sample_sizes

#model with symmetric migration from Moments bitbucket
def split_mig_moments(params, ns, pop_ids=None):
    if pop_ids is not None and len(pop_ids) != 2:
        raise ValueError("pop_ids must be a list of two population IDs")
    nu1, nu2, T, m, theta1 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1],
                                                  theta= theta1)
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T, m=np.array([[0, m], [m, 0]]), theta= theta1)
    fs.pop_ids = pop_ids
    return fs

func_moments = split_mig_moments

#boundaries for nu1, nu2, Ts, m12, theta1
upper_bound = [10, 10, 5, 50, 5e5]
lower_bound = [1e-5, 1e-5, 1e-5, 0, 3.75e5]

#constants
mu = 2.8e-9 #from Keightley et al. 2014
L = int(L_file) #imported from system argument #2
g = 0.07692308 #equals 13 gen/year. Calculated based on biological intuition.
theta_mean = 419800 #mean theta when multiplying 0.005*L, where Watterson's theta=0.005
theta_sd = 41980 #sd for theta when multiplying 0.0005*L, where 0.0005 is the standard deviation of Watterson's theta

# run X optimizations from uniform-sampled starting params
for i in range(int(iterations)): #iterations is imported from sys. argument #1
    print("starting optimization "+str(i))

    #Start the run by picking random parameters from a uniform distribution for nu1, nu2, Ts, and m
    #theta has own distribution from which it is pulled (Gaussian), where the mean is the: (average Watterson's pi)*L
    popt=[np.random.uniform(lower_bound[0], upper_bound[0]), np.random.uniform(lower_bound[1], upper_bound[1]), np.random.uniform(lower_bound[2], upper_bound[2]), np.random.uniform(lower_bound[3], upper_bound[3]), np.random.normal(loc=theta_mean, scale=theta_sd, size=None)]

    #This is the optimization step for moments.
    #popt is the prior.
    #fs folded is a tranform SFS by folding it. The original SFS loaded in sys.arg #1 is quasi-folded. i.e. polarized to reference genome.
    #Folding it is a must, because reference is not 100% ancestral
    #multinom is set to false as theta is being estimated in the model as a parameter.
    popt=moments.Inference.optimize_log(popt, fs_folded, func_moments,
                                        lower_bound=lower_bound, upper_bound=upper_bound,
                                        verbose=False, maxiter=100, multinom=False,
                                        )

    #This number is 5. i.e., count parameters. there is an opportunity to streamline the code by propagating this from the beginning.
    params = len(["nu1", "nu2", "Ts", "m12", "theta1"]) #for use in AIC calculation

    #This is the moments function.
    model = func_moments(popt, ns)

    #Calculate log likelihood of the model fit
    #ll is calculated using Poisson approach instead of calculating model theta (as it is a parameter here)
    #from the moments manual (6.1 pg. 17): "In the Poisson approach, the likelihood is the product of Poisson likelihoods for each 
    #entry in the data FS, given an expected value from the model FS"
    ll_model= moments.Inference.ll(model, fs_folded)
    #Now calculate AIC of model fit
    aic = 2*params - 2*ll_model
    print('Maximum log composite likelihood: {0}'.format(ll_model))

    #Now estimate theta from model fit
    #As theta is a scalar for the SFS, and we're already parameterizing this scalar, this should equal ~1
    theta_model = moments.Inference.optimal_sfs_scaling(model, fs_folded)
    
    #reducing complexity of calculations to follow by adding variables in lieu of expressions/esoteric list calls
    nu1=popt[0]
    nu2=popt[1]
    Ts=popt[2]
    m12=popt[3]
    theta_param=popt[4]
    Nref= theta_param/(4*mu*L)

    #Now calculate Ts from Model fit
    divergence_time = 2*Nref*Ts*g #calculates divergence time in years

    #Now calculate Migration rate (fraction of migrants that move between pops)
    Mij = m12/(2*Nref) #actual migration rate
    #Below is an old code which had an error. Keric has since updated it. kept for record keeping. #Jcbn Jun18,2021
    #mig_pop1 = Mij*(2*popt[0]) #number of individuals going i to j
    #mig_pop2 = Mij*(2*popt[1]) #number of individuals going j to i

    #Now we are estimated the nominal migration rate based on Mij
    mig_pop1 = Mij*(nu1*Nref) #number of individuals going i to j: migrants = Mij*nu1*Nref
    mig_pop2 = Mij*(nu2*Nref) #number of individuals going j to i: migrants = Mij*nu2*Nref

    #Now estimate population size
    pop1_size = nu1*Nref #pop1 size
    pop2_size = nu2*Nref #pop2 size

    #Open the output file
    PMmod=open('./outputs_param/%s_output.txt' % Pair_name,'a')

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
        str(theta_model)+'\t'+ #theta as calculated by the model- sanity check-- should equal ~1
        str(theta_param)+'\t'+ #theta as calculated as a parameter. should equal theta model, just doing sanity check
        str(nu1)+'\t'+ #raw parameter output
        str(nu2)+'\t'+ #raw parameter output
        str(Ts)+'\t'+ #raw parameter output
        str(m12)+'\t'+ #raw parameter output
        str(S)+'\t'+ #gives number of segregating sites used in sfs... sanity check
        str(ll_model)+'\t'+
        str(aic)+'\n')
    PMmod.close()

print("Moments finished running")