#import packages that'll be used
import moments
from moments import Numerics
from moments import Integration
from moments import Spectrum
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

#define sys args
#DEBUGGED AS FOR LOOPS
fs_file = sys.argv[1]
L_file = sys.argv[2]
iterations = sys.argv[3]
Pair_name = sys.argv[4]

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
data = moments.Spectrum.from_file(fs_file)
fs_folded = data.fold()
ns = fs_folded.sample_sizes

# These are the grid point settings will use for extrapolation.
# They're somewhat arbitrary, but should be big enough that they don't change LL if increased
pts_l = [100,140,180]

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
L = int(L_file) #specific to Austria dataset
g = 0.0383562 #equals 2 weeks converted to years

# run X optimizations from uniform-sampled starting params
for i in range(int(iterations)):
    print("starting optimization "+str(i))
    popt=[np.random.uniform(lower_bound[x],upper_bound[x]) for x in range(4)]
    popt=moments.Inference.optimize_log(popt, fs_folded, func_moments,
                                        lower_bound=lower_bound, upper_bound=upper_bound,
                                        verbose=False, maxiter=100,
                                        )
    params = len(["nu1", "nu2", "Ts", "m"]) #for use in AIC calculation
    model = func_moments(popt, ns)
    ll_model=moments.Inference.ll_multinom(model, fs_folded)
    aic = 2*params - 2*ll_model
    print('Maximum log composite likelihood: {0}'.format(ll_model))
    theta = moments.Inference.optimal_sfs_scaling(model, fs_folded)
    divergence_time = 2*(theta/(4*mu*L))*popt[2]*g #calculates divergence time in years
    Mij = popt[3]/(2*(theta/(4*mu*L))) #actual migration rate
    mig_pop1 = Mij*(2*popt[0]) #number of individuals going i to j
    mig_pop2 = Mij*(2*popt[1]) #number of individuals going j to i
    pop1_size = popt[0]*(theta/(4*mu*L)) #pop1 size
    pop2_size = popt[1]*(theta/(4*mu*L)) #pop2 size

    PMmod=open('%s_output.txt' % Pair_name,'a')

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
