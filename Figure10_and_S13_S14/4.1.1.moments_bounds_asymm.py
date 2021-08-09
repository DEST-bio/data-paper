# This script will run the program moments by Jouganous, J., Long, W., Ragsdale, A. P., & Gravel, S. (2017)
# Parameters are constrained at 1e-5 for this run 
# Written by Keric Lamb, UVA 2021
# ksl2za@virginia.edu
# import packages that'll be used
import moments
from moments import Numerics
from moments import Integration
from moments import Spectrum
import dadi
from dadi import Misc
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
from datetime import datetime

#define sys args
#DEBUGGED AS FOR LOOPS
fs_file = sys.argv[1]
L_file = sys.argv[2]
iterations = sys.argv[3]
Pair_name = sys.argv[4]
pop_name1 = sys.argv[6]
pop_name2 = sys.argv[5]
pool_n1 = sys.argv[7]
pool_n2 = sys.argv[8]

#converting floats to integers
pool_n1= int(pool_n1)
pool_n2= int(pool_n2)

#read in data as file
dd = dadi.Misc.make_data_dict(fs_file) #reads in genomalicious SNP file
data = pd.read_csv(fs_file, sep="\t", nrows=1)

#setting up if else check 
#if names fed by metadat match dd then will run as expected
#else if names are not equal, it swaps them
#leaves pool_n1 alone
if pop_name1==data.columns[3]:
    pop_id1=pop_name1
    pop_id2=pop_name2
else:
    pop_id1=pop_name2
    pop_id2=pop_name1

#setting pop id's and projections from if/else
pop_id=[pop_id1,pop_id2]
pools=[pool_n1, pool_n2]

#makes sfs
fs_folded = Spectrum.from_data_dict(dd, pop_ids=pop_id, projections=pools, polarized=False) #takes data dict and folds
ns = fs_folded.sample_sizes #gets sample sizes of dataset
S = fs_folded.S()
#fs_folded.mask[:1,:] = True
#fs_folded.mask[ :,:1] = True

#gets time for warning file
now = datetime.now()
#generates warning file
PMmod=open("../warnings.txt", 'a')
PMmod.write(
    str("%s" % now)+'\t'+ 
    str("%s" % Pair_name)+'\t'+
    str("bound")+'\t'+
    str("%s" % S)+'\n')
PMmod.close()

#tells if sfs is empty and if so, kills job
if S==0:
    quit()
else:
    print("continuing")

#opening output file to give column names
PMmod=open('%s_output.bound.asymm.txt' % Pair_name,'w')
PMmod.write(
            str("Pair_name")+'\t'+ #print pair name
            str("fs_name")+'\t'+ #double checking fs_lines[y] is working as I want it to
            str("L")+'\t'+ #double checking L is working as I want it to
            str("pop1_size")+'\t'+ #nu1
            str("pop2_size")+'\t'+ #nu2
            str("divergence_time")+'\t'+ #divergence T
            str("mig_pop1")+'\t'+ #Migration ij
            str("mig_pop2")+'\t'+ #Migration ji
            str("Mij")+'\t'+
            str("Mji")+'\t'+
            str("theta")+'\t'+
            str("nu1")+'\t'+
            str("nu2")+'\t'+
            str("Ts")+'\t'+
            str("m12")+'\t'+
            str("m21")+'\t'+
            str("fs_sanitycheck")+'\t'+
            str("-2LL_model")+'\t'+
            str("AIC")+'\n')
PMmod.close()

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
    nu1, nu2, T, m12, m21 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.integrate([nu1, nu2], T, m=np.array([[0, m12], [m21, 0]]))
    fs.pop_ids = pop_ids
    return fs

func_moments = split_mig_moments

#boundaries
#bounded to same standards as theta param script
upper_bound = [10, 10, 5, 50, 50]
lower_bound = [1e-5, 1e-5, 1e-5, 0, 0]

#constants
mu = 2.8e-9 #from Keightley et al. 2014
L = int(L_file) #imported from system argument #2
g = 0.07692308 #equals 13 gen/year. Calculated based on biological intuition.

# run X optimizations from uniform-sampled starting params
for i in range(int(iterations)): #iterations is imported from sys. argument #1
    print("starting optimization "+str(i))

    #Parameters are set above " nu1, nu2, T, m12, m21"
    #This number is 5. i.e., count parameters. 
    params = len(["nu1", "nu2", "Ts", "m12", "m21"]) #for use in AIC calculation
    
    #Start the run by picking random parameters from a uniform distribution.
    popt=[np.random.uniform(lower_bound[x],upper_bound[x]) for x in range(params)]

    #This is the optimization step for moments.
    #popt is the prior.
    #fs folded is a tranform SFS by folding it. The original SFS loaded in sys.arg #1 is quasi-folded. i.e. polarized to reference genome.
    #Folding it is a must, because reference is not 100% ancestral
    popt=moments.Inference.optimize_log(popt, fs_folded, func_moments,
                                        lower_bound=lower_bound, upper_bound=upper_bound,
                                        verbose=False, maxiter=100,
                                        )

    #This is the moments function.
    model = func_moments(popt, ns)

    #Calculate log likelihood of the model fit
    ll_model=moments.Inference.ll_multinom(model, fs_folded)
    #Now calculate AIC of model fit
    aic = 2*params - 2*ll_model
    print('Maximum log composite likelihood: {0}'.format(ll_model))

    #Now estimate theta from model fit
    theta = moments.Inference.optimal_sfs_scaling(model, fs_folded)
    
    #reducing complexity of calculations to follow by adding variables in lieu of expressions/esoteric df calls
    Nref= theta/(4*mu*L)
    nu1=popt[0]
    nu2=popt[1]
    Ts=popt[2]
    m12=popt[3]
    m21=popt[4]
    
    #Now calculate Ts from Model fit
    divergence_time = 2*Nref*Ts*g #calculates divergence time in years

    #Now calculate Migration rate (fraction of migrants that move between pops)
    Mij = m12/(2*Nref) #actual migration rate for populations from 2 into 1
    Mji = m21/(2*Nref) #actual migration rate for populations from 1 into 2
    
    #Now we are estimated the nominal migrant count based on Mij
    mig_pop1 = Mij*(nu1*Nref) #number of individuals going i to j: migrants = Mij*nu1*Nref
    mig_pop2 = Mji*(nu2*Nref) #number of individuals going j to i: migrants = Mij*nu2*Nref

    #Now estimate population size
    pop1_size = nu1*Nref #pop1 size
    pop2_size = nu2*Nref #pop2 size

    #Open the output file
    PMmod=open('%s_output.bound.asymm.txt' % Pair_name,'a')

    #Dumping output ot outfile
    PMmod.write(
        str(Pair_name)+'\t'+ #print pair name
        str(fs_file)+'\t'+ #double checking fs is the right one
        str(L)+'\t'+ #double checking L is working as desired
        str(pop1_size)+'\t'+ #nu1
        str(pop2_size)+'\t'+ #nu2
        str(divergence_time)+'\t'+ #divergence T
        str(mig_pop1)+'\t'+ #Migrants ij
        str(mig_pop2)+'\t'+ #Migrants ji
        str(Mij)+'\t'+ #migration rate ij
        str(Mji)+'\t'+ #migration rate ji
        str(theta)+'\t'+
        str(nu1)+'\t'+
        str(nu2)+'\t'+
        str(Ts)+'\t'+
        str(m12)+'\t'+
        str(m21)+'\t'+
        str(S)+'\t'+ #sanity check... should give number of segregating sites in SFS
        str(ll_model)+'\t'+
        str(aic)+'\n')
    PMmod.close()

print("Moments finished running")