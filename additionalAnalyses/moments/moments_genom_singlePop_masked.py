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
import pandas as pd

#define sys args
#DEBUGGED AS FOR LOOPS
fs_file = sys.argv[1]
L_file = sys.argv[2]
Pair_name = sys.argv[3]
pop_id1 = sys.argv[4]
pop_id2 = sys.argv[5]
projection1 = sys.argv[6]
projection2 = sys.argv[7]

#need to make projection numbers into integers. Python doesn't read in sys arg's as ints.
projection1= int(projection1)
projection2= int(projection2)

#opening output file to give column names
PMmod=open('%s_onepop_output.txt' % Pair_name,'w')
PMmod.write(
            str("Pair_name")+'\t'+ #print pair name
            str("fs_name")+'\t'+ #double checking fs_lines[y] is working as I want it to
            str("L")+'\t'+ #double checking L is working as I want it to
    		str("Nref")+'\t'+ #Nref pop size
            str("theta")+'\t'+
            str("-2LL_model")+'\t'+
            str("AIC")+'\n')
PMmod.close()

#read in data as file
dd = Misc.make_data_dict(fs_file) #reads in genomalicious SNP file

#making pop_id and projection lists for making the folded spectrum
pop_id=[pop_id1,pop_id2]
projection=[projection1,projection2]

fs_folded = Spectrum.from_data_dict(dd, pop_ids=pop_id, projections=projection, polarized=False) #takes data dict and folds
ns = fs_folded.sample_sizes #gets sample sizes of dataset
fs_folded.mask[:1,:] = True
fs_folded.mask[ :,:1] = True

#this code was used to import Alan's custom fs files. Now commented out as we've switched to genomalicious (R package) compiled SNP files
# data = moments.Spectrum.from_file(fs_file)
# fs_folded = data.fold()
# ns = fs_folded.sample_sizes

#standard neutral model from Moments 2D demography manual example
def snm(ns, pop_ids=pop_id):
    """
    ns = [n1, n2]

    Standard neutral model with a split but no divergence.

    :param ns: List of population sizes in first and second populations.
    :param pop_ids: List of population IDs.
    """
    if pop_ids is not None and len(pop_ids) != 2:
        raise ValueError("pop_ids must be a list of two population IDs")
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])
    fs.pop_ids = pop_ids
    return fs

func_moments = snm

# run model
print("starting optimization")

params = 0 #for use in AIC calculation. 0 as there are no parameters for the snm

#This is the moments function.
model = func_moments(ns)

mu = 2.8e-9 #from Keightley et al. 2014
L = int(L_file) #imported from system argument #2
g = 0.07692308 #equals 13 gen/year. Calculated based on biological intuition.

ll_model=moments.Inference.ll_multinom(model, fs_folded)
#Now calculate AIC of model fit
aic = 2*params - 2*ll_model
print('Maximum log composite likelihood: {0}'.format(ll_model))

#Now estimate theta from model fit
theta = moments.Inference.optimal_sfs_scaling(model, fs_folded)
#Now estimate population size
Nref = theta/(4*mu*L) #pop size

#Open the output file
PMmod=open('%s_onepop_output.txt' % Pair_name,'a')

#Dumping output ot outfile
PMmod.write(
    str(Pair_name)+'\t'+ #print pair name
    str(fs_file)+'\t'+ #double checking fs is the right one
    str(L)+'\t'+ #double checking L is working as desired
    str(Nref)+'\t'+ #Nref pop size
    str(theta)+'\t'+
    str(ll_model)+'\t'+
    str(aic)+'\n')
PMmod.close()

print("Moments finished running")
