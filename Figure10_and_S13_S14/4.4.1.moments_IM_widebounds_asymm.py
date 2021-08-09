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
print("done reading in dd and check")

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

#generate sfs
fs_folded = Spectrum.from_data_dict(dd, pop_ids=pop_id, projections=pools, polarized=False) #takes data dict and folds
print("done converting dd to fs")
ns = fs_folded.sample_sizes #gets sample sizes of dataset
S = fs_folded.S()
#fs_folded.mask[:1,:] = True
#fs_folded.mask[ :,:1] = True

#generate time
now = datetime.now()
#generates warning file
PMmod=open("../warnings.txt", 'a')
PMmod.write(
    str("%s" % now)+'\t'+ 
    str("%s" % Pair_name)+'\t'+
    str("theta")+'\t'+
    str("%s" % S)+'\n')
PMmod.close()

#if file is not working (sfs is empty), this kills the script
if S==0:
    quit()
else:
    print("continuing")

#opening output file to give column names
PMmod=open('%s_output.IM.widebounds.asymm.txt' % Pair_name,'w')
PMmod.write(
    str("Pair_name")+'\t'+ #print pair name
    str("fs_file")+'\t'+ #double checking fs is the right one
    str("L")+'\t'+ #double checking L is working as desired
    str("pop1_bottle_size")+'\t'+ #size of pop1 after split and instantaneous size change
    str("pop2_bottle_size")+'\t'+ #size of pop2 after split and instantaneous size change
    str("pop1_final_size")+'\t'+ #nu1 converted to inds
    str("pop2_final_size")+'\t'+ #nu2 converted to inds
    str("divergence_time")+'\t'+ #divergence T
    str("Mij")+'\t'+ #Migration ij
    str("Mji")+'\t'+
    str("mig_pop1")+'\t'+ #Migrants in pop1 from pop2
    str("mig_pop2")+'\t'+ #Migrants in pop2 from pop1
    str("theta")+'\t'+ #theta as calculated by the model- sanity check-- should equal ~1
    str("nu1B")+'\t'+ #raw parameter output
    str("nu1F")+'\t'+ #raw parameter output
    str("nu2B")+'\t'+ #raw parameter output
    str("nu2F")+'\t'+ #raw parameter output
    str("Ts")+'\t'+ #raw parameter output
    str("m12")+'\t'+ #raw parameter output
    str("m21")+'\t'+
    str("S")+'\t'+ #gives number of segregating sites used in sfs... sanity check
    str("ll_model")+'\t'+
    str("aic")+'\n')
PMmod.close()

#model with asymmetric growth migration from Moments bitbucket with 1 change: instantaneous change rather than
#pop sizes post split of s and s-1 (now s1 and s2)
def IM_bottle_growth(params, ns):
    """
    Model with split, bottleneck, exp recovery, migration
    nu1B: The bottleneck size for pop1
    nu2F: The final size for pop1
    nu2B: The bottleneck size for pop2
    nu2F: The final size for pop2
    m: The scaled migration rate
    T: The time between the split and present
    theta1: theta parameter for the model to optimize

    ns = n1,n2: Size of fs to generate.
    """
    nu1B, nu1F, nu2B, nu2F, T, m12, m21 = params

    # f for the equilibrium ancestral population
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])

    # We need to define a function to describe the non-constant population 2
    # size. lambda is a convenient way to do so.
    nu1_func = lambda t: nu1B * (nu1F / nu1B) ** (t / T)
    nu2_func = lambda t: nu2B * (nu2F / nu2B) ** (t / T)
    nu_func = lambda t: [nu1_func(t), nu2_func(t)]
    fs.integrate(nu_func, T, m=np.array([[0, m12], [m21, 0]]))
    
    return fs

func_moments = IM_bottle_growth

#boundaries for nu1B, nu1F, nu2B, nu2F, Ts, m12, m21
upper_bound = [10, 10, 10, 10, 5, 50, 50]
lower_bound = [1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 0, 0]

#constants
mu = 2.8e-9 #from Keightley et al. 2014
L = int(L_file) #imported from system argument #2
g = 0.07692308 #equals 13 gen/year. Calculated based on biological intuition.

# run X optimizations from uniform-sampled starting params
for i in range(int(iterations)): #iterations is imported from sys. argument #1
    print("starting optimization "+str(i))

    #This number is 7. i.e., count parameters.
    params = len(["nu1B", "nu1F", "nu2B", "nu2F", "Ts", "m12", "m21"]) #for use in AIC calculation

    #Start the run by picking random parameters from a uniform distribution for nu1, nu2, Ts, and m12, m21
    #theta has own distribution from which it is pulled (Gaussian), where the mean is the: (average Watterson's pi)*L
    popt=[np.random.uniform(lower_bound[x],upper_bound[x]) for x in range(params)]

    #This is the optimization step for moments.
    #popt is the prior.
    #fs folded is a tranform SFS by folding it. The original SFS loaded in sys.arg #1 is quasi-folded. i.e. polarized to reference genome.
    #Folding it is a must, because reference is not 100% ancestral
    #multinom is set to false as theta is being estimated in the model as a parameter.
    popt=moments.Inference.optimize_log(popt, fs_folded, func_moments,
                                        lower_bound=lower_bound, upper_bound=upper_bound,
                                        verbose=True, maxiter=3,
                                        )

    #This is the moments function.
    model = func_moments(popt, ns)

    #Calculate log likelihood of the model fit
    #ll is calculated using Poisson approach instead of calculating model theta (as it is a parameter here)
    #from the moments manual (6.1 pg. 17): "In the Poisson approach, the likelihood is the product of Poisson likelihoods for each 
    #entry in the data FS, given an expected value from the model FS"
    ll_model= moments.Inference.ll_multinom(model, fs_folded)
    #Now calculate AIC of model fit
    aic = 2*params - 2*ll_model
    print('Maximum log composite likelihood: {0}'.format(ll_model))

    #Now estimate theta from model fit
    #As theta is a scalar for the SFS, and we're already parameterizing this scalar, this should equal ~1
    theta = moments.Inference.optimal_sfs_scaling(model, fs_folded)
    
    #reducing complexity of calculations to follow by adding variables in lieu of expressions/esoteric list calls
    nu1B=popt[0]
    nu1F=popt[1]
    nu2B=popt[2]
    nu2F=popt[3]
    Ts=popt[4]
    m12=popt[5]
    m21=popt[6]
    Nref= theta/(4*mu*L)

    #Now calculate Ts from Model fit
    divergence_time = 2*Nref*Ts*g #calculates divergence time in years

    #Now calculate Migration rate (fraction of migrants that move between pops)
    Mij = m12/(2*Nref) #actual migration rate from pop2 into pop1
    Mji = m21/(2*Nref) #actual migration rate from pop2 into pop1

    #Now we are estimated the nominal migrant count based on Mij
    mig_pop1 = Mij*(nu1F*Nref) #number of individuals in pop1 from pop2
    mig_pop2 = Mji*(nu2F*Nref) #number of individuals in pop2 from pop1

    #Now estimate population size
    pop1_bottle_size = nu1B*Nref #size of pop1 after split and instantaneous size change
    pop2_bottle_size = nu2B*Nref #size of pop2 after split and instantaneous size change
    pop1_final_size = nu1F*Nref #pop1 size
    pop2_final_size = nu2F*Nref #pop2 size

    #Open the output file
    PMmod=open('%s_output.IM.widebounds.asymm.txt' % Pair_name,'a')

    #Dumping output ot outfile
    PMmod.write(
        str(Pair_name)+'\t'+ #print pair name
        str(fs_file)+'\t'+ #double checking fs is the right one
        str(L)+'\t'+ #double checking L is working as desired
        str(pop1_bottle_size)+'\t'+ #size of pop1 after split and instantaneous size change
        str(pop2_bottle_size)+'\t'+ #size of pop2 after split and instantaneous size change
        str(pop1_final_size)+'\t'+ #nu1 converted to inds
        str(pop2_final_size)+'\t'+ #nu2 converted to inds
        str(divergence_time)+'\t'+ #divergence T
        str(Mij)+'\t'+ #Migration ij
        str(Mji)+'\t'+#Migration ji
        str(mig_pop1)+'\t'+ #Migrants in pop1 from pop2
        str(mig_pop2)+'\t'+ #Migrants in pop2 from pop1
        str(theta)+'\t'+ #theta as calculated by the model
        str(nu1B)+'\t'+ #raw parameter output
        str(nu1F)+'\t'+ #raw parameter output
        str(nu2B)+'\t'+ #raw parameter output
        str(nu2F)+'\t'+ #raw parameter output
        str(Ts)+'\t'+ #raw parameter output
        str(m12)+'\t'+ #raw parameter output
        str(m21)+'\t'+ #raw parameter output
        str(S)+'\t'+ #gives number of segregating sites used in sfs... sanity check
        str(ll_model)+'\t'+
        str(aic)+'\n')
    PMmod.close()

print("Moments finished running")