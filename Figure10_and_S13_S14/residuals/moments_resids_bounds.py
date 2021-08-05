#script to print residual graphs
#import libraries
#works for any run where theta is not parameterized
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
import  pandas as pd
import pylab
import datetime
from datetime import datetime

#setting system argv's
fs_file = sys.argv[1]
Pair_name = sys.argv[2]
pop_name1 = sys.argv[3]
pop_name2 = sys.argv[4]
pool_n1 = sys.argv[5]
pool_n2 = sys.argv[6]
nu1 = float(sys.argv[7])
nu2 = float(sys.argv[8])
m12 = float(sys.argv[9])
Ts = float(sys.argv[10])
theta = float(sys.argv[11])
resid = int(sys.argv[12])

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

fs_folded = Spectrum.from_data_dict(dd, pop_ids=pop_id, projections=pools, polarized=False) #takes data dict and folds
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

# Calculate the best-fit model AFS.
popt= [nu1, nu2, Ts, m12]
model = func_moments(popt, ns)
# Likelihood of the data given the model AFS.
ll_model = moments.Inference.ll_multinom(model, fs_folded)
print("Maximum log composite likelihood: {0}".format(ll_model))
# The optimal value of theta given the model.
theta = moments.Inference.optimal_sfs_scaling(model, fs_folded)
print("Optimal value of theta: {0}".format(theta))
    
fig = pylab.figure(1)
fig.clear()
moments.Plotting.plot_2d_comp_multinom(
    model, fs_folded, resid_range=resid, pop_ids=('%s' % pop_id1, '%s' % pop_id2), adjust=True
)
fig.savefig("%s_bounds_resids.png" % Pair_name, dpi=100)