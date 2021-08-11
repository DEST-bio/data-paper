# This script will run the program moments by Jouganous, J., Long, W., Ragsdale, A. P., & Gravel, S. (2017)
# Script to output residuals from metadata file
# Generalized so it can run any model we ran, provided proper metadata input
# Symmetric models script can process: theta split_mig, theta IMbg, wide bounds split_mig, wide bounds IMbg
# Asymmetric models script can process: wide bounds split_mig, wide bounds IMbg
# Written by Keric Lamb, UVA 2021
# ksl2za@virginia.edu

#import relevant packages
import moments
from moments import Numerics
from moments import Integration
from moments import Spectrum
from moments import Inference
import dadi
from dadi import Misc #make data dict step was not working in moments for unknown reason
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')  # need this when DISPLAY is not defined (i.e. slurm runs)
import matplotlib.pyplot as plt
import  pandas as pd
import pylab
import datetime
from datetime import datetime

#setting system argv's
fs_file = sys.argv[1]
Pair_name = sys.argv[2]
param = sys.argv[3]
model_type= sys.argv[4]
model_sym= sys.argv[5]
pop_name1 = sys.argv[6]
pop_name2 = sys.argv[7]
pool_n1 = sys.argv[8]
pool_n2 = sys.argv[9]

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

#if file is not working (sfs is empty), this kills the run
if S==0:
    quit()
else:
    print("continuing")

#begin model and param residual runs

if param==str("theta"):
    if model_type==str("split_mig"):
        #bring in system inputs
        nu1 = float(sys.argv[10])
        nu2 = float(sys.argv[11])
        m12 = float(sys.argv[12])
        Ts = float(sys.argv[14])
        theta = float(sys.argv[15]) #note this should be theta_param, not theta_model

        #generate time
        now = datetime.now()
        #generates warning file
        PMmod=open("warnings.txt", 'a')
        PMmod.write(
            str("%s" % now)+'\t'+
            str("%s" % Pair_name)+'\t'+
            str("symmetric")+'\t'+
            str("split_mig")+'\t'+
            str("theta")+'\t'+
            str("%s" % S)+'\n')
        PMmod.close()

        #model with symmetric migration from Moments bitbucket
        #theta is parameterized here
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

        func_moments = split_mig_moments #generalizes function name for ease of c/p

        # Calculate the best-fit model AFS.
        popt= [nu1, nu2, Ts, m12, theta]
        model = func_moments(popt, ns)
        # Likelihood of the data given the model AFS.
        ll_model= moments.Inference.ll(model, fs_folded)
        print("Maximum log composite likelihood: {0}".format(ll_model))

        fig = pylab.figure(1)
        fig.clear()
        moments.Plotting.plot_2d_comp_Poisson(
            model, fs_folded, resid_range=10, pop_ids=('%s' % pop_id1, '%s' % pop_id2), adjust=True
        )
        fig.savefig("%s_theta_splitmig_resids.png" % Pair_name, dpi=100)

        #grab residual array from internal moments function
        resid = moments.Inference.Anscombe_Poisson_residual(model, fs_folded)
        resid = pd.DataFrame(resid) #convert the spectrum residual object to a df
        resid.to_csv("%s_theta_splitmig_resids_array.csv" % Pair_name, sep='\t', index = True)

    else:
        #bring in system inputs
        m12 = float(sys.argv[12])
        Ts = float(sys.argv[14])
        theta = float(sys.argv[15]) #note this should be theta_param, not theta_model
        nu1B = float(sys.argv[16])
        nu2B = float(sys.argv[17])
        nu1F = float(sys.argv[18])
        nu2F = float(sys.argv[19])

        #generate time
        now = datetime.now()
        #generates warning file
        PMmod=open("warnings.txt", 'a')
        PMmod.write(
            str("%s" % now)+'\t'+
            str("%s" % Pair_name)+'\t'+
            str("symmetric")+'\t'+
            str("IMbg")+'\t'+
            str("theta")+'\t'+
            str("%s" % S)+'\n')
        PMmod.close()

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
            nu1B, nu1F, nu2B, nu2F, T, m, theta1 = params

            # f for the equilibrium ancestral population
            sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1], theta=theta1)
            fs = moments.Spectrum(sts)
            fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])

            # We need to define a function to describe the non-constant population 2
            # size. lambda is a convenient way to do so.
            nu1_func = lambda t: nu1B * (nu1F / nu1B) ** (t / T)
            nu2_func = lambda t: nu2B * (nu2F / nu2B) ** (t / T)
            nu_func = lambda t: [nu1_func(t), nu2_func(t)]
            fs.integrate(nu_func, T, m=np.array([[0, m], [m, 0]]), theta=theta1)

            return fs

        func_moments = IM_bottle_growth

        # Calculate the best-fit model AFS.
        popt= [nu1B, nu1F, nu2B, nu2F, Ts, m12, theta]
        model = func_moments(popt, ns)
        # Likelihood of the data given the model AFS.
        ll_model= moments.Inference.ll(model, fs_folded)
        print("Maximum log composite likelihood: {0}".format(ll_model))

        fig = pylab.figure(1)
        fig.clear()
        moments.Plotting.plot_2d_comp_Poisson(
            model, fs_folded, resid_range=10, pop_ids=('%s' % pop_id1, '%s' % pop_id2), adjust=True
        )
        fig.savefig("%s_theta_IMbg_resids.png" % Pair_name, dpi=100)

        #grab residual array from internal moments function
        resid = moments.Inference.Anscombe_Poisson_residual(model, fs_folded)
        resid = pd.DataFrame(resid) #convert the spectrum residual object to a df
        resid.to_csv("%s_theta_IMbg_resids_array.csv" % Pair_name, sep='\t', index = True)

else: #non-theta models
    if model_sym==str("symm"):
        if model_type==str("split_mig"):

            #import parameters relevant to split_mig widebound models
            nu1 = float(sys.argv[10])
            nu2 = float(sys.argv[11])
            m12 = float(sys.argv[12])
            Ts = float(sys.argv[14])

            #generate time
            now = datetime.now()
            #generates warning file
            PMmod=open("warnings.txt", 'a')
            PMmod.write(
                str("%s" % now)+'\t'+
                str("%s" % Pair_name)+'\t'+
                str("symmetric")+'\t'+
                str("split_mig")+'\t'+
                str("widebounds")+'\t'+
                str("%s" % S)+'\n')
            PMmod.close()

            #model with symmetric migration from Moments bitbucket
            #theta is parameterized here
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
            ll_model= moments.Inference.ll_multinom(model, fs_folded)
            print("Maximum log composite likelihood: {0}".format(ll_model))

            fig = pylab.figure(1)
            fig.clear()
            moments.Plotting.plot_2d_comp_multinom(
                model, fs_folded, resid_range=10, pop_ids=('%s' % pop_id1, '%s' % pop_id2), adjust=True
            )
            fig.savefig("%s_widebounds_splitmig_resids.png" % Pair_name, dpi=100)

            #grab residual array from internal moments function
            resid = moments.Inference.Anscombe_Poisson_residual(model, fs_folded)
            resid = pd.DataFrame(resid) #convert the spectrum residual object to a df
            resid.to_csv("%s_widebounds_splitmig_resids_array.csv" % Pair_name, sep='\t', index = True)

        else: #symmetric widebounds IMbg
            #bring in system inputs
            m12 = float(sys.argv[12])
            Ts = float(sys.argv[14])
            nu1B = float(sys.argv[16])
            nu2B = float(sys.argv[17])
            nu1F = float(sys.argv[18])
            nu2F = float(sys.argv[19])

            #generate time
            now = datetime.now()
            #generates warning file
            PMmod=open("warnings.txt", 'a')
            PMmod.write(
                str("%s" % now)+'\t'+
                str("%s" % Pair_name)+'\t'+
                str("symmetric")+'\t'+
                str("IMbg")+'\t'+
                str("widebounds")+'\t'+
                str("%s" % S)+'\n')
            PMmod.close()

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
                nu1B, nu1F, nu2B, nu2F, T, m = params

                # f for the equilibrium ancestral population
                sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
                fs = moments.Spectrum(sts)
                fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1])

                # We need to define a function to describe the non-constant population 2
                # size. lambda is a convenient way to do so.
                nu1_func = lambda t: nu1B * (nu1F / nu1B) ** (t / T)
                nu2_func = lambda t: nu2B * (nu2F / nu2B) ** (t / T)
                nu_func = lambda t: [nu1_func(t), nu2_func(t)]
                fs.integrate(nu_func, T, m=np.array([[0, m], [m, 0]]))

                return fs

            func_moments = IM_bottle_growth

            # Calculate the best-fit model AFS.
            popt= [nu1B, nu1F, nu2B, nu2F, Ts, m12]
            model = func_moments(popt, ns)
            # Likelihood of the data given the model AFS.
            ll_model= moments.Inference.ll_multinom(model, fs_folded)
            print("Maximum log composite likelihood: {0}".format(ll_model))

            fig = pylab.figure(1)
            fig.clear()
            moments.Plotting.plot_2d_comp_multinom(
                model, fs_folded, resid_range=10, pop_ids=('%s' % pop_id1, '%s' % pop_id2), adjust=True
            )
            fig.savefig("%s_widebounds_IMbg_resids.png" % Pair_name, dpi=100)

            #grab residual array from internal moments function
            resid = moments.Inference.Anscombe_Poisson_residual(model, fs_folded)
            resid = pd.DataFrame(resid) #convert the spectrum residual object to a df
            resid.to_csv("%s_widebounds_IMbg_resids_array.csv" % Pair_name, sep='\t', index = True)

    else: #asymmetric runs
        if model_type==str("split_mig"):
            #bring in system inputs
            nu1 = float(sys.argv[10])
            nu2 = float(sys.argv[11])
            m12 = float(sys.argv[12])
            m21 = float(sys.argv[13])
            Ts = float(sys.argv[14])

            #generate time
            now = datetime.now()
            #generates warning file
            PMmod=open("warnings.txt", 'a')
            PMmod.write(
                str("%s" % now)+'\t'+
                str("%s" % Pair_name)+'\t'+
                str("asymmetric")+'\t'+
                str("split_mig")+'\t'+
                str("widebounds")+'\t'+
                str("%s" % S)+'\n')
            PMmod.close()

            #model with asymmetric migration from Moments bitbucket
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

            # Calculate the best-fit model AFS.
            popt= [nu1, nu2, Ts, m12, m21]
            model = func_moments(popt, ns)
            # Likelihood of the data given the model AFS.
            ll_model= moments.Inference.ll_multinom(model, fs_folded)
            print("Maximum log composite likelihood: {0}".format(ll_model))

            fig = pylab.figure(1)
            fig.clear()
            moments.Plotting.plot_2d_comp_multinom(
                model, fs_folded, resid_range=10, pop_ids=('%s' % pop_id1, '%s' % pop_id2), adjust=True
            )
            fig.savefig("%s_widebounds_splitmig_asymm_resids.png" % Pair_name, dpi=100)

            #grab residual array from internal moments function
            resid = moments.Inference.Anscombe_Poisson_residual(model, fs_folded)
            resid = pd.DataFrame(resid) #convert the spectrum residual object to a df
            resid.to_csv("%s_widebounds_splitmig_asymm_resids_array.csv" % Pair_name, sep='\t', index = True)

        else: #asymmetric IMbg widebounds model
            #bring in system inputs
            m12 = float(sys.argv[12])
            m21 = float(sys.argv[13])
            Ts = float(sys.argv[14])
            nu1B = float(sys.argv[16])
            nu2B = float(sys.argv[17])
            nu1F = float(sys.argv[18])
            nu2F = float(sys.argv[19])

            #generate time
            now = datetime.now()
            #generates warning file
            PMmod=open("warnings.txt", 'a')
            PMmod.write(
                str("%s" % now)+'\t'+
                str("%s" % Pair_name)+'\t'+
                str("asymmetric")+'\t'+
                str("IMbg")+'\t'+
                str("widebounds")+'\t'+
                str("%s" % S)+'\n')
            PMmod.close()

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
                nu1B, nu1F, nu2B, nu2F, T, m12, m21= params

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

            # Calculate the best-fit model AFS.
            popt= [nu1B, nu1F, nu2B, nu2F, Ts, m12, m21]
            model = func_moments(popt, ns)
            # Likelihood of the data given the model AFS.
            ll_model= moments.Inference.ll_multinom(model, fs_folded)
            print("Maximum log composite likelihood: {0}".format(ll_model))

            #plot the figure
            fig = pylab.figure(1)
            fig.clear()
            moments.Plotting.plot_2d_comp_multinom(
                model, fs_folded, resid_range=10, pop_ids=('%s' % pop_id1, '%s' % pop_id2), adjust=True
            )
            fig.savefig("%s_widebounds_IMbg_asymm_resids.png" % Pair_name, dpi=100)

            #grab residual array from internal moments function


            resid = moments.Inference.Anscombe_Poisson_residual(model, fs_folded)
            resid = pd.DataFrame(resid) #convert the spectrum residual object to a df
            resid.to_csv("%s_widebounds_IMbg_asymm_resids_array.csv" % Pair_name, sep='\t', index = True)

print("residual script is done running %s" % Pair_name)
