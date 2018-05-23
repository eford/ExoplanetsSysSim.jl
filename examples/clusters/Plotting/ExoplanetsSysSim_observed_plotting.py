#To import required modules:
import numpy as np
import time
import matplotlib
import matplotlib.cm as cm #for color maps
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec #for specifying plot attributes
from matplotlib import ticker #for setting contour plots to log scale
import scipy.integrate #for numerical integration
import scipy.misc #for factorial function
#matplotlib.rc('text', usetex=True)





##### This module will be used to plot simulated data generated from the ExoplanetsSysSim (Eric's version of the clustering method)

#To define some useful constants:
N_Kep = 150061 #number of Kepler targets satisfying our cuts to give our observed catalog

AU = 1.496*10.**13. #AU in cm
Msun = 1.989*10.**30. #Solar mass in kg
Rsun = 6.957*10.**10. #Solar radius in cm
Rearth = 6.371*10.**8. #Earth radius in cm

savefigures = False
loadfiles_directory = 'Clustering_Method_Figures/ExoplanetsSysSim/Non_clustered/Sim4/' #'ExoplanetsSysSim.jl-master/examples/clusters/'; 'Clustering_Method_Figures/ExoplanetsSysSim/Power_law_r1_r2_sigma_r/'
savefigures_directory = 'Clustering_Method_Figures/ExoplanetsSysSim/Non_clustered/'
run_number = ''
model_name = 'ExoplanetsSysSim_Non_clustered_Model4' + run_number #'ExoplanetsSysSim_Non_clustered_Model'

param_keys_all = [("num_targets_sim_pass_one", r'$N_{\rm stars,sim}$'),
                  ("log_rate_clusters", r'$\lambda_c$'),
                  ("max_clusters_in_sys", r'$N_{c,\rm max}$'),
                  #("log_rate_planets_per_cluster", r'$\lambda_p$'),
                  #("max_planets_in_clusters", r'$N_{p,\rm max}$'),
                  ("power_law_P", r'$\alpha_P$'),
                  ("min_period", r'$P_{\rm min}$'),
                  ("max_period", r'$P_{\rm max}$'),
                  #("power_law_r", r'$\alpha_R$'),
                  ("power_law_r1", r'$\alpha_{R1}$'),
                  ("power_law_r2", r'$\alpha_{R2}$'),
                  ("min_radius (R_earth)", r'$R_{p,\rm min}$ $(R_\oplus)$'),
                  ("max_radius (R_earth)", r'$R_{p,\rm max}$ $(R_\oplus)$'),
                  ("break_radius (R_earth)", r'$R_{p,\rm break}$ $(R_\oplus)$'),
                  ("sigma_incl", r'$\sigma_i$ (deg)'),
                  ("sigma_incl_near_mmr", r'$\sigma_{i,\rm res}$ (deg)'),
                  ("sigma_hk", r'$\sigma_e$'),
                  ("num_mutual_hill_radii", r'$\Delta_c$'),
                  ("mr_power_index", r'$\alpha_{mr}$'),
                  ("mr_max_mass (M_earth)", r'$M_{p,\rm max}$ $(M_\oplus)$'),
                  #("sigma_log_radius_in_cluster", r'$\sigma_R$'),
                  #("sigma_logperiod_per_pl_in_cluster", r'$\sigma_N$')
                  ] #list of the symbols and names for all the model parameters; NOTE: although the params are named log rate of clusters and planets per cluster, we use the symbols and values for the rates





##### To load the files with the systems with observed planets:

#To first read the number of simulated targets and bounds for the periods and radii:
with open(loadfiles_directory + 'periods%s.out' % run_number, 'r') as file:
    for line in file:
        if line[:26] == '# num_targets_sim_pass_one':
            N_sim = int(line[28:])
        elif line[:14] == '# max_incl_sys':
            max_incl_sys = float(line[16:])
            cos_factor = np.cos(max_incl_sys*np.pi/180.)
        elif line[:12] == '# min_period':
            P_min = float(line[14:])
        elif line[:12] == '# max_period':
            P_max = float(line[14:])
        elif line[:12] == '# min_radius':
            radii_min = float(line[24:])
        elif line[:12] == '# max_radius':
            radii_max = float(line[24:])

#To read the simulation parameters from the file:
param_vals_all = [] #list to be filled with the values of all the model parameters
with open(loadfiles_directory + 'periods%s.out' % run_number, 'r') as file:
    for line in file:
        for i in range(len(param_keys_all)):
            chars = len(param_keys_all[i][0])
            if line[:3+chars] == '# ' + param_keys_all[i][0] + ':':
                if param_keys_all[i][0][:3] == 'log':
                    param_vals_all.append(np.round(np.exp(float(line[4+chars:])), 4))
                elif (param_keys_all[i][0][:11] == 'num_targets') or (param_keys_all[i][0][:11] == 'mr_max_mass'):
                    param_vals_all.append(int(float(line[4+chars:])))
                else:
                    param_vals_all.append(np.round(float(line[4+chars:]), 4))
if len(param_vals_all) != len(param_keys_all):
    print 'Problem with reading parameter values...'



P_per_sys = [] #list to be filled with lists of the observed periods per system (days)
with open(loadfiles_directory + 'periods%s.out' % run_number, 'r') as file:
    for line in file:
        if line[0] != '#':
            line = line[1:-2]
            line_per_sys = line.split('; ')
            print len(line_per_sys)
            for x in line_per_sys:
                P_sys = x.split()
                P_sys = [float(i) for i in P_sys]
                P_per_sys.append(P_sys)
                #print P_sys

D_per_sys = [] #list to be filled with lists of the transit depths per system
with open(loadfiles_directory + 'depths%s.out' % run_number, 'r') as file:
    for line in file:
        if line[0] != '#':
            line = line[1:-2]
            line_per_sys = line.split('; ')
            #print len(line_per_sys)
            for x in line_per_sys:
                D_sys = x.split()
                D_sys = [float(i) for i in D_sys]
                D_per_sys.append(D_sys)

tdur_per_sys = [] #list to be filled with lists of the transit durations per system (days)
with open(loadfiles_directory + 'durations%s.out' % run_number, 'r') as file:
    for line in file:
        if line[0] != '#':
            line = line[1:-2]
            line_per_sys = line.split('; ')
            #print len(line_per_sys)
            for x in line_per_sys:
                tdur_sys = x.split()
                tdur_sys = [float(i) for i in tdur_sys]
                tdur_per_sys.append(tdur_sys)

Mstar_obs = [] #list to be filled with the stellar masses of the systems with observed planets (Msun)
with open(loadfiles_directory + 'stellar_masses_obs%s.out' % run_number, 'r') as file:
    for line in file:
        if line[0] != '#':
            line = line[1:-2]
            Mstars = line.split(', ')
            Mstars = [float(i) for i in Mstars]
            Mstar_obs = Mstar_obs + Mstars
Mstar_obs = np.array(Mstar_obs)

Rstar_obs = [] #list to be filled with the stellar radii of the systems with observed planets (Rsun)
with open(loadfiles_directory + 'stellar_radii_obs%s.out' % run_number, 'r') as file:
    for line in file:
        if line[0] != '#':
            line = line[1:-2]
            Rstars = line.split(', ')
            Rstars = [float(i) for i in Rstars]
            Rstar_obs = Rstar_obs + Rstars
Rstar_obs = np.array(Rstar_obs)

P_obs = [] #list to be zero-padded so each list of periods is sorted and has the same length, and then converted to an array
D_obs = [] #list to be zero-padded so each list of depths is sorted (by period) and has the same length, and then converted to an array
tdur_obs = [] #list to be zero-padded so each list of transit durations is sorted (by period) and has the same length, and then converted to an array

Pmin = 0. #set a minimum period (days), discarding planets less than this period

Mmax = len(P_per_sys[-1]) #maximum planet multiplicity generated by the clustering method
for i in range(len(P_per_sys)):
    i_sorted = np.argsort(P_per_sys[i]) #array of indices which would sort the system by period
    P_sorted = np.array(P_per_sys[i])[i_sorted]
    P_sorted_cut = P_sorted[P_sorted > Pmin]
    D_sorted_cut = np.array(D_per_sys[i])[i_sorted][P_sorted > Pmin]
    tdur_sorted_cut = np.array(tdur_per_sys[i])[i_sorted][P_sorted > Pmin]
    
    P_sys = list(P_sorted_cut) + [0]*(Mmax - len(P_sorted_cut)) #zero-pad the list up to Mmax elements
    D_sys = list(D_sorted_cut) + [0]*(Mmax - len(D_sorted_cut)) #zero-pad the list up to Mmax elements
    tdur_sys = list(tdur_sorted_cut) + [0]*(Mmax - len(tdur_sorted_cut)) #zero-pad the list up to Mmax elements

    P_obs.append(P_sys)
    D_obs.append(D_sys)
    tdur_obs.append(tdur_sys)
P_obs = np.array(P_obs)
D_obs = np.array(D_obs)
tdur_obs = np.array(tdur_obs)*24.*60. #tdur_obs converted to mins

Mtot_obs = np.sum(P_obs > 0, axis=1) #array of observed planet multiplicites
Nmult_obs = np.array([np.sum(Mtot_obs == x) for x in range(1,Mmax+1)]) #array of total numbers of systems with observed planet multiplicities of 1,2,3,...,Mmax planets
radii_obs = np.sqrt(D_obs)*np.transpose([Rstar_obs])*(Rsun/Rearth) #array of planet radii, in Earth radii



#To calculate the observed period ratios, period-normalized transit duration ratios, and transit depth ratios:
res_ratios, res_width = [1.5, 2.0], 0.05 #NOTE: in the model, the near-resonant planets have period ratios between X and (1+w)*X where X = [2/1, 3/2, 4/3, 5/4] and w = 0.05!

Rm_obs = [] #list to be filled with the observed period ratios
D_ratio_obs = [] #list to be filled with the observed transit depth ratios
xi_obs = [] #list to be filled with the period-normalized transit duration ratios
xi_res_obs = [] #list to be filled with the period-normalized transit duration ratios for planet pairs near resonance
xi_res32_obs = []
xi_res21_obs = []
xi_nonres_obs = [] #list to be filled with the period_normalized transit duration ratios for planet pairs not in resonance
for i in range(len(P_obs)):
    P_obs_system = P_obs[i][P_obs[i] > 0]
    tdur_obs_system = tdur_obs[i][P_obs[i] > 0]
    D_obs_system = D_obs[i][P_obs[i] > 0]
    
    #To calculate all the observed period ratios:
    Rm_obs_system = list(P_obs_system[1:]/P_obs_system[0:-1]) #list of period ratios observed in this system
    Rm_obs_system = np.array(Rm_obs_system + [0]*(Mmax - 1 - len(Rm_obs_system))) #to add filler 0's to Rm_obs_system to pad it to Mmax - 1 elements
    Rm_obs.append(Rm_obs_system)
    
    #To calculate all the observed transit depth ratios:
    D_ratio_obs_system = list(D_obs_system[1:]/D_obs_system[0:-1]) #list of transit depth ratios observed in this system
    D_ratio_obs_system = np.array(D_ratio_obs_system + [0]*(Mmax - 1 - len(D_ratio_obs_system))) #to add filler 0's to D_ratio_obs_system to pad it to Mmax - 1 elements
    D_ratio_obs.append(D_ratio_obs_system)
    
    #To calculate all the period-normalized transit duration ratios:
    xi_obs_system = list((tdur_obs_system[0:-1]/tdur_obs_system[1:])*(P_obs_system[1:]/P_obs_system[0:-1])**(1./3.)) #list of period-normalized transit duration ratios in this system
    xi_obs_system = np.array(xi_obs_system + [0]*(Mmax - 1 - len(xi_obs_system))) #to add filler 0's to xi_obs_system to pad it to Mmax - 1 elements
    xi_obs.append(xi_obs_system)
    
    #To separate the period-normalized transit duration ratios for planet pairs near vs. not in resonance:
    mask_res_system = np.zeros(len(Rm_obs_system), dtype=bool)
    mask_res32_system = np.zeros(len(Rm_obs_system), dtype=bool)
    mask_res21_system = np.zeros(len(Rm_obs_system), dtype=bool)
    
    mask_res_system[(Rm_obs_system >= res_ratios[0]) & (Rm_obs_system <= res_ratios[0]*(1.+res_width))] = 1
    mask_res_system[(Rm_obs_system >= res_ratios[1]) & (Rm_obs_system <= res_ratios[1]*(1.+res_width))] = 1
    mask_res32_system[(Rm_obs_system >= res_ratios[0]) & (Rm_obs_system <= res_ratios[0]*(1.+res_width))] = 1
    mask_res21_system[(Rm_obs_system >= res_ratios[1]) & (Rm_obs_system <= res_ratios[1]*(1.+res_width))] = 1
    
    xi_res_obs_system = list(xi_obs_system[mask_res_system])
    xi_res32_obs_system = list(xi_obs_system[mask_res32_system])
    xi_res21_obs_system = list(xi_obs_system[mask_res21_system])
    xi_nonres_obs_system = list(xi_obs_system[~mask_res_system])
    xi_res_obs_system = np.array(xi_res_obs_system + [0]*(10 - len(xi_res_obs_system)))
    xi_res32_obs_system = np.array(xi_res32_obs_system + [0]*(10 - len(xi_res32_obs_system)))
    xi_res21_obs_system = np.array(xi_res21_obs_system + [0]*(10 - len(xi_res21_obs_system)))
    xi_nonres_obs_system = np.array(xi_nonres_obs_system + [0]*(10 - len(xi_nonres_obs_system)))
    xi_res_obs.append(xi_res_obs_system)
    xi_res32_obs.append(xi_res32_obs_system)
    xi_res21_obs.append(xi_res21_obs_system)
    xi_nonres_obs.append(xi_nonres_obs_system)

Rm_obs = np.array(Rm_obs)
D_ratio_obs = np.array(D_ratio_obs)
xi_obs = np.array(xi_obs)
xi_res_obs = np.array(xi_res_obs)
xi_res32_obs = np.array(xi_res32_obs)
xi_res21_obs = np.array(xi_res21_obs)
xi_nonres_obs = np.array(xi_nonres_obs)

P_obs_flat = P_obs.flatten() #all the observed periods of all the planets
P_obs_flat = P_obs_flat[P_obs_flat > 0]

Rm_obs_flat = Rm_obs.flatten() #all the observed period ratios of all the observed adjacent planets
Rm_obs_flat = Rm_obs_flat[Rm_obs_flat > 0]

D_obs_flat = D_obs.flatten() #all the transit depths
D_obs_flat = D_obs_flat[D_obs_flat > 0]
radii_obs_flat = radii_obs.flatten() #all the observed planet radii, in Earth radii
radii_obs_flat = radii_obs_flat[radii_obs_flat > 0]

D_ratio_obs_flat = D_ratio_obs.flatten() #all the transit depth ratios
D_ratio_obs_flat = D_ratio_obs_flat[D_ratio_obs_flat > 0]

tdur_obs_flat = tdur_obs.flatten() #all the observed transit durations, in mins
tdur_obs_flat = tdur_obs_flat[tdur_obs_flat > 0]

xi_obs_flat = xi_obs.flatten() #all the observed period-normalized transit duration ratios
xi_obs_flat = xi_obs_flat[xi_obs_flat > 0]
xi_res_obs_flat = xi_res_obs.flatten() #the observed period-normalized transit duration ratios for planet pairs near resonance
xi_res_obs_flat = xi_res_obs_flat[xi_res_obs_flat > 0]
xi_res32_obs_flat = xi_res32_obs.flatten()
xi_res32_obs_flat = xi_res32_obs_flat[xi_res32_obs_flat > 0]
xi_res21_obs_flat = xi_res21_obs.flatten()
xi_res21_obs_flat = xi_res21_obs_flat[xi_res21_obs_flat > 0]
xi_nonres_obs_flat = xi_nonres_obs.flatten() #the observed period-normalized transit duration ratios for planet pairs not in resonance
xi_nonres_obs_flat = xi_nonres_obs_flat[xi_nonres_obs_flat > 0]





#To plot the resulting 'observed' distributions (observed total planet multiplicity, periods, period ratios):
'''
fig = plt.figure(figsize=(16,8))
plot = GridSpec(2,2,left=0.075,bottom=0.115,right=0.95,top=0.925,wspace=0.15,hspace=0.4)

ax = plt.subplot(plot[0,0])
plt.title('Observed planet multiplicity', fontsize=24)
x = Mtot_obs[Mtot_obs > 0]
counts, bins = np.histogram(x, bins=np.max(x)+1, range=(-0.5, np.max(x)+0.5))
bins_mid = (bins[:-1] + bins[1:])/2.
plt.plot(bins_mid, counts, 'o-', color='k', label='%s systems with observed planets' % len(x))
ax.tick_params(axis='both', labelsize=20)
plt.xlim([np.min(x), np.max(x)])
plt.xlabel(r'$M_{\rm tot}$', fontsize=20)
plt.ylabel('Number', fontsize=20)
plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[0,1])
plt.title('Observed periods', fontsize=24)
hist = plt.hist(P_obs_flat, bins=np.logspace(np.log10(np.min(P_obs_flat)), np.log10(np.max(P_obs_flat)), 21), histtype='step', log=True, color='k', label=r'Observed planets')
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([np.min(P_obs_flat)-1, 1.1*np.max(P_obs_flat)])
plt.xlabel(r'$P$ (days)', fontsize=20)
#plt.ylabel('Number', fontsize=20)
#plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[1,0])
plt.title('Observed period ratios ($P_{i+1}/P_i < 5$)', fontsize=24)
x = Rm_obs_flat[Rm_obs_flat < 5]
plt.hist(x, bins=np.linspace(1,5,51), histtype='step', color='k', label='Observed adjacent pairs')
#plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([1,5])
plt.xlabel(r'$P_{i+1}/P_i$', fontsize=20)
plt.ylabel('Number', fontsize=20)
#plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[1,1])
plt.title('Observed period ratios ($P_{i+1}/P_i > 5$)', fontsize=24)
x = Rm_obs_flat[Rm_obs_flat > 5]
plt.hist(x, bins=50, histtype='step', color='k', label='Observed adjacent pairs')
#plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([5,np.max(x)])
plt.xlabel(r'$P_{i+1}/P_i$', fontsize=20)
#plt.ylabel('Number', fontsize=20)
#plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

if savefigures == True:
    plt.savefig(savefigures_directory + model_name + '_observed.pdf')
else:
    plt.show()
plt.close()
'''


#To plot the resulting 'observed' distributions (observed transit durations, period-normalized transit duration ratios, transit depths, and transit depth ratios):
'''
fig = plt.figure(figsize=(16,8))
plot = GridSpec(2,2,left=0.075,bottom=0.115,right=0.95,top=0.925,wspace=0.15,hspace=0.4)

ax = plt.subplot(plot[0,0])
plt.title('Transit durations', fontsize=24)
x = tdur_obs_flat
#x = x[x < 1000.]
plt.hist(x, bins=50, histtype='step', color='k')
ax.tick_params(axis='both', labelsize=20)
plt.xlim([np.min(x), np.max(x)])
plt.xlabel(r'$t_{\rm dur}$ (mins)', fontsize=20)
plt.ylabel('Number', fontsize=20)

ax = plt.subplot(plot[0,1])
plt.title('Period-normalized transit duration ratios', fontsize=24)
x = np.log10(xi_obs_flat)
hist = plt.hist(x, bins=50, histtype='step', color='k', label='All')
ax.tick_params(axis='both', labelsize=20)
plt.xlim([np.min(x), np.max(x)])
plt.xlabel(r'$\log{\/xi}$', fontsize=20)
#plt.ylabel('Number', fontsize=20)
#plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[1,0])
plt.title('Transit depths', fontsize=24)
x = D_obs_flat
plt.hist(x, bins=np.logspace(np.log10(np.min(x)), np.log10(np.max(x)), 51), histtype='step', color='k')
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([np.min(x), np.max(x)])
plt.xlabel(r'$\delta$', fontsize=20)
plt.ylabel('Number', fontsize=20)
#plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[1,1])
plt.title('Transit depth ratios', fontsize=24)
x = D_ratio_obs_flat
hist = plt.hist(x, bins=np.logspace(np.log10(np.min(x)), np.log10(np.max(x)), 51), histtype='step', color='k')
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([np.min(x), np.max(x)])
plt.xlabel(r'$\delta_{i+1}/\delta_i$', fontsize=20)
#plt.ylabel('Number', fontsize=20)
#plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

if savefigures == True:
    plt.savefig(savefigures_directory + model_name + '_observed_properties.pdf')
else:
    plt.show()
plt.close()
'''


#To plot the resulting 'observed' distributions (observed period-normalized transit duration ratios):
'''
fig = plt.figure(figsize=(16,8))
plot = GridSpec(2,2,left=0.075,bottom=0.115,right=0.95,top=0.925,wspace=0.15,hspace=0.4)

ax = plt.subplot(plot[0,0])
plt.title('Period-normalized transit duration ratios', fontsize=24)
x = np.log10(xi_obs_flat)
hist = plt.hist(x, bins=50, histtype='step', color='k', label='All')
plt.hist(np.log10(xi_res_obs_flat), bins=hist[1], histtype='step', color='m', label='Near resonance')
plt.hist(np.log10(xi_nonres_obs_flat), bins=hist[1], histtype='step', color='g', label='Not near resonance')
ax.tick_params(axis='both', labelsize=20)
plt.xlim([np.min(x), np.max(x)])
plt.xlabel(r'$\log{\/xi}$', fontsize=20)
plt.ylabel('Number', fontsize=20)
plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[1,0])
plt.title('Period-normalized transit duration ratios', fontsize=24)
x = np.log10(xi_obs_flat)
hist = plt.hist(x, bins=np.linspace(-0.5,0.5,51), histtype='step', weights=np.ones(len(x))/len(x), color='k', label='All')
plt.hist(np.log10(xi_res_obs_flat), bins=hist[1], histtype='step', weights=np.ones(len(xi_res_obs_flat))/float(len(xi_res_obs_flat)), color='m', label='Near resonance')
plt.hist(np.log10(xi_nonres_obs_flat), bins=hist[1], histtype='step', weights=np.ones(len(xi_nonres_obs_flat))/float(len(xi_nonres_obs_flat)), color='g', label='Not near resonance')
ax.tick_params(axis='both', labelsize=20)
plt.xlim([hist[1][0], hist[1][-1]])
plt.xlabel(r'$\log{\/xi}$', fontsize=20)
plt.ylabel('Fraction', fontsize=20)
#plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[1,1])
plt.title('Period-normalized transit duration ratios', fontsize=24)
plt.hist(np.log10(xi_res32_obs_flat), bins=hist[1], histtype='step', weights=np.ones(len(xi_res32_obs_flat))/float(len(xi_res32_obs_flat)), color='r', label='Near 3:2 MMR')
plt.hist(np.log10(xi_res21_obs_flat), bins=hist[1], histtype='step', weights=np.ones(len(xi_res21_obs_flat))/float(len(xi_res21_obs_flat)), color='b', label='Near 2:1 MMR')
ax.tick_params(axis='both', labelsize=20)
plt.xlim([hist[1][0], hist[1][-1]])
plt.xlabel(r'$\log{\/xi}$', fontsize=20)
#plt.ylabel('Fraction', fontsize=20)
plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

if savefigures == True:
    plt.savefig(savefigures_directory + model_name + '_observed_logxi.pdf')
else:
    plt.show()
plt.close()
'''





##### To compare the simulated observed distributions to the Kepler observed distributions using the K-S distance:

#To define functions for computing distances between two distributions:
def CRPD_dist(En, On): #NOTE: this distance can potentially give negative values?!
    #This function computes the Cressie Read Power Divergence statistic for observed planet multiplicities
    #En and On must be arrays of the total numbers of systems with 1,2,3,... observed planets, in the simulated (i.e. expected) and the actual (i.e. observed Kepler) data, respectively
    E_array = En/float(np.sum(En)) #normalized numbers (fractions) of simulated systems with 1,2,3,... observed planets
    O_array = On/float(np.sum(On)) #normalized numbers (fractions) of actual Kepler systems with 1,2,3,... observed planets
    rho = 0.
    for i,E_i in enumerate(E_array):
        if En[i] != 0:
            rho += O_array[i]*((O_array[i]/E_array[i])**(2./3.) - 1.)
    rho = (9./5.)*rho
    return rho

def KS_dist_mult(x1, x2):
    #This function computes the K-S distance for two discrete, integer distributions (for multiplicities)
    #This function returns two values: the K-S distance and the x value corresponding to that distance
    x12_max = np.max((np.max(x1), np.max(x2))) #maximum value of x1 and x2
    x1_counts, x1_bins = np.histogram(x1, bins=x12_max, range=(0.5, x12_max+0.5))
    x2_counts, x2_bins = np.histogram(x2, bins=x12_max, range=(0.5, x12_max+0.5))
    pdf_diffs = x1_counts/np.float(len(x1)) - x2_counts/np.float(len(x2))
    cdf_diffs = np.cumsum(pdf_diffs)
    KS_dist = np.max(np.abs(cdf_diffs)) #K-S distance
    KS_x = np.arange(1, x12_max+1)[np.where(np.abs(cdf_diffs) == KS_dist)[0][0]] #x value where the distance is the largest
    return KS_dist, KS_x

def KS_dist(x1, x2):
    #This function computes the K-S distance for two continuous distributions (no repeated values)
    #This function returns two values: the K-S distance and the x value corresponding to that distance
    x_all = np.concatenate((x1, x2)) #combined array
    i_all_sorted = np.argsort(x_all) #array of indices that would sort the combined array
    pdf_diffs = np.concatenate((np.ones(len(x1))/np.float(len(x1)), -np.ones(len(x2))/np.float(len(x2))))[i_all_sorted]
    cdf_diffs = np.cumsum(pdf_diffs)
    KS_dist = np.max(np.abs(cdf_diffs)) #K-S distance
    KS_x = x_all[i_all_sorted][np.where(np.abs(cdf_diffs) == KS_dist)[0][0]] #x value (a value in either x1 or x2) where the distance is the largest
    return KS_dist, KS_x

def AD_dist(x1, x2):
    #This function computes and returns the AD distance for two continuous distributions (no repeated values), according to A. N. Pettitt (1976) Eq. (1.2)
    n, m = len(x1), len(x2)
    N = n + m
    x_all = np.concatenate((x1, x2)) #combined array
    i_all_sorted = np.argsort(x_all) #array of indices that would sort the combined array
    M_i_diffs = np.concatenate((np.ones(n), np.zeros(m)))[i_all_sorted]
    M_i_array = np.cumsum(M_i_diffs)[:-1] #array of M_i except for last element, i.e. from i=1 to i=N-1
    i_array = 1. + np.arange(N-1) #array of i from i=1 to i=N-1
    AD_dist = (1./(n*m))*np.sum(((M_i_array*N - n*i_array)**2.)/(i_array*(N - i_array))) #AD distance
    return AD_dist

def AD_dist2(x1, x2): #I tested this and it returns the same results as AD_dist()
    #This function computes and returns the AD distance for two continuous distributions (no repeated values), according to Scholz & Stephens (1987) Eq. (3)
    n1, n2 = len(x1), len(x2)
    N = n1 + n2
    x_all = np.concatenate((x1, x2)) #combined array
    i_all_sorted = np.argsort(x_all) #array of indices that would sort the combined array
    
    M_1j_diffs = np.concatenate((np.ones(n1), np.zeros(n2)))[i_all_sorted]
    M_1j_array = np.cumsum(M_1j_diffs)[:-1] #array of M_1j except for last element, i.e. from j=1 to j=N-1
    M_2j_diffs = np.concatenate((np.zeros(n1), np.ones(n2)))[i_all_sorted]
    M_2j_array = np.cumsum(M_2j_diffs)[:-1] #array of M_2j except for last element, i.e. from j=1 to j=N-1
    j_array = 1. + np.arange(N-1) #array of j from j=1 to j=N-1
    
    AD_dist = (1./N)*((1./n1)*np.sum(((N*M_1j_array - n1*j_array)**2.)/(j_array*(N - j_array))) + (1./n2)*np.sum(((N*M_2j_array - n2*j_array)**2.)/(j_array*(N - j_array)))) #AD distance
    return AD_dist


#To load and compute the exoplanet multiplicities, periods, and period ratios of the confirmed Kepler exoplanets:
Q1Q17_DR25 = np.genfromtxt('q1_q17_dr25_koi.tab_selectcols_new.csv', dtype={'names': ('KepID', 'KOI', 'Archive_Disp', 'Kepler_Disp', 'Disp', 'P', 't_D', 'depth', 'Rp', 'Rstar'), 'formats': ('i8', 'S9', 'S15', 'S15', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8',)}, delimiter=',', usecols=(1,2,3,4,5,10,11,12,13,14)) #orbit periods 'P' are in days; transit durations 't_D' are in hrs; transit depths 'depth' are in ppm; planetary radii 'Rp' are in Rearth; stellar radii 'Rstar' are in Rsolar
Q1Q17_DR25 = Q1Q17_DR25[1:] #skip_header doesn't work so manually get rid of first row of NaNs

Q1Q17_DR25_stellar = np.genfromtxt('q1_q17_dr25_stellar_koi.tab_selectcols.csv', dtype={'names': ('KepID', 'mag', 'teff', 'logg', 'cdpp1_5', 'cdpp2', 'cdpp2_5', 'cdpp3', 'cdpp3_5', 'cdpp4_5', 'cdpp5', 'cdpp6', 'cdpp7_5', 'cdpp9', 'cdpp10_5', 'cdpp12', 'cdpp12_5', 'cdpp15'), 'formats': ('i8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8',)}, delimiter=',', usecols=(2,3,4,5,12,13,14,15,16,17,18,19,20,21,22,23,24,25)) #RMS CDPP's are in ppm?
Q1Q17_DR25_stellar = Q1Q17_DR25_stellar[1:] #skip_header doesn't work so manually get rid of first row of NaNs

table_confirmed = Q1Q17_DR25[(Q1Q17_DR25['Archive_Disp'] == 'CONFIRMED') + (Q1Q17_DR25['Archive_Disp'] == 'CANDIDATE')]
table_stellar = Q1Q17_DR25_stellar

#To make cuts in period and planetary radii:
table_confirmed = table_confirmed[(table_confirmed['P'] > P_min) & (table_confirmed['P'] < P_max)]
table_confirmed = table_confirmed[(table_confirmed['Rp'] > radii_min) & (table_confirmed['Rp'] < radii_max)]

#To make cuts based on stellar properties of T_eff and logg (and CDPP if we choose to):
teff_confirmed = np.zeros(len(table_confirmed)) #list to be filled with the T_eff (K) for the objects
logg_confirmed = np.zeros(len(table_confirmed)) #list to be filled with the logg(cgs) for the objects
cdpp5_confirmed = np.zeros(len(table_confirmed)) #list to be filled with the RMS CDPP 5h values for the objects
cdpp6_confirmed = np.zeros(len(table_confirmed)) #list to be filled with the RMS CDPP 6h values for the objects
for i,KepID in enumerate(table_confirmed['KepID']):
    teff_confirmed[i] = table_stellar['teff'][table_stellar['KepID'] == KepID]
    logg_confirmed[i] = table_stellar['logg'][table_stellar['KepID'] == KepID]
    cdpp5_confirmed[i] = table_stellar['cdpp5'][table_stellar['KepID'] == KepID]
    cdpp6_confirmed[i] = table_stellar['cdpp6'][table_stellar['KepID'] == KepID]

cdpp_cut = 250.
print 'Fraction of CONFIRMED and CANDIDATE planets after CDPP cut: %s/%s' % (len(table_confirmed[(teff_confirmed > 4000.) & (teff_confirmed < 7000.) & (logg_confirmed > 4.) & (cdpp5_confirmed < cdpp_cut)]), len(table_confirmed[(teff_confirmed > 4000.) & (teff_confirmed < 7000.) & (logg_confirmed > 4.)]))

table_confirmed = table_confirmed[(teff_confirmed > 4000.) & (teff_confirmed < 7000.) & (logg_confirmed > 4.) & (cdpp5_confirmed < cdpp_cut)]


#To compute the planet multiplicities and period ratios:
KOI_systems = np.array([x[:6] for x in table_confirmed['KOI']])
checked_bools = np.zeros(len(table_confirmed)) #0's denote KOI that were not checked yet; 1's denote already checked KOI
M_confirmed = [] #list to be filled with the planet multiplicities of systems with confirmed planets
R_confirmed = [] #list to be filled with period ratios of adjacent confirmed planet pairs

P_sys_confirmed = [] #list to be filled with lists of planet periods per system
radii_sys_confirmed = [] #list to be filled with lists of planet radii per system

D_ratio_confirmed = [] #list to be filled with the transit depth ratios of adjacent confirmed planet pairs
xi_confirmed = [] #list to be filled with the period-normalized transit duration ratios of adjacent confirmed planet pairs
xi_res_confirmed = [] #list to be filled with the period-normalized transit duration ratios of adjacent confirmed planet pairs near resonance
xi_res32_confirmed = [] #list to be filled with the period-normalized transit duration ratios of adjacent confirmed planet pairs near 3:2 resonance
xi_res21_confirmed = [] #list to be filled with the period-normalized transit duration ratios of adjacent confirmed planet pairs near 2:1 resonance
xi_nonres_confirmed = [] #list to be filled with the period-normalized transit duration ratios of adjacent confirmed planet pairs not in resonance
t_D_confirmed = table_confirmed['t_D'] #array of the transit durations (hrs) of all the confirmed planets
D_confirmed = table_confirmed['depth']/(1e6) #array of the transit depths (fraction) of all the confirmed planets
radii_confirmed = table_confirmed['Rp'] #array of the planetary radii (Rearth) of all the confirmed planets
D_confirmed = D_confirmed[radii_confirmed < 10.]
for i in range(len(KOI_systems)):
    if checked_bools[i] == 0: #if the KOI has not been checked (included while looking at another planet in the same system)
        system_i = np.where(KOI_systems == KOI_systems[i])[0]
        checked_bools[system_i] = 1
        
        #To get the periods and transit durations in this system:
        system_P = table_confirmed['P'][system_i] #periods of all the planets in this system
        system_t_D = table_confirmed['t_D'][system_i] #transit durations of all the planets in this system
        system_D = table_confirmed['depth'][system_i] #transit depths of all the planets in this system
        system_radii = table_confirmed['Rp'][system_i] #radii of all the planets in this system
        system_sort_i = np.argsort(system_P) #indices that would sort the periods of the planets in this system
        system_P = system_P[system_sort_i] #periods of all the planets in this system, sorted
        system_t_D = system_t_D[system_sort_i] #transit durations of all the planets in this system, sorted by period
        system_D = system_D[system_sort_i] #transit depths of all the planets in this system, sorted by period
        system_radii = system_radii[system_sort_i] #radii of all the planets in this system, sorted by period
        
        #To count the total number of planets in this system:
        M_confirmed.append(len(system_P))
        P_sys_confirmed.append(list(system_P) + [0]*(6 - len(system_P)))
        radii_sys_confirmed.append(list(system_radii) + [0]*(6 - len(system_radii)))
        
        #To compute the period ratios and period-normalized transit duration ratios in this system (and separate into planet pairs near vs. not in resonance):
        system_R = system_P[1:]/system_P[0:-1] #period ratios of all the adjacent planet pairs in this system
        system_D_ratio = system_D[1:]/system_D[0:-1] #transit depth ratios of all the adjacent planet pairs in this system
        system_xi = (system_t_D[0:-1]/system_t_D[1:])*(system_P[1:]/system_P[0:-1])**(1./3.) #period-normalized transit duration ratios of all the adjacent planet pairs in this system
        mask_res_system = np.zeros(len(system_R), dtype=bool)
        mask_res32_system = np.zeros(len(system_R), dtype=bool)
        mask_res21_system = np.zeros(len(system_R), dtype=bool)
        
        mask_res_system[(system_R >= res_ratios[0]) & (system_R <= res_ratios[0]*(1.+res_width))] = 1
        mask_res_system[(system_R >= res_ratios[1]) & (system_R <= res_ratios[1]*(1.+res_width))] = 1
        mask_res32_system[(system_R >= res_ratios[0]) & (system_R <= res_ratios[0]*(1.+res_width))] = 1
        mask_res21_system[(system_R >= res_ratios[1]) & (system_R <= res_ratios[1]*(1.+res_width))] = 1
        system_xi_res = system_xi[mask_res_system]
        system_xi_res32 = system_xi[mask_res32_system]
        system_xi_res21 = system_xi[mask_res21_system]
        system_xi_nonres = system_xi[~mask_res_system]
        #if sum(mask_res_system) > 0:
        #print system_R[mask_res_system], system_xi_res
        
        for R in system_R:
            R_confirmed.append(R)
        for D_ratio in system_D_ratio:
            D_ratio_confirmed.append(D_ratio)
        for xi in system_xi:
            xi_confirmed.append(xi)
        for xi in system_xi_res:
            xi_res_confirmed.append(xi)
        for xi in system_xi_res32:
            xi_res32_confirmed.append(xi)
        for xi in system_xi_res21:
            xi_res21_confirmed.append(xi)
        for xi in system_xi_nonres:
            xi_nonres_confirmed.append(xi)

P_sys_confirmed = np.array(P_sys_confirmed)
radii_sys_confirmed = np.array(radii_sys_confirmed)

P_confirmed = table_confirmed['P']
M_confirmed = np.array(M_confirmed)
Nmult_confirmed = np.array([np.sum(M_confirmed == x) for x in range(1,np.max((Mmax,np.max(M_confirmed)))+1)])
R_confirmed = np.array(R_confirmed)
D_ratio_confirmed = np.array(D_ratio_confirmed)
xi_confirmed = np.array(xi_confirmed)
xi_res_confirmed = np.array(xi_res_confirmed)
xi_res32_confirmed = np.array(xi_res32_confirmed)
xi_res21_confirmed = np.array(xi_res21_confirmed)
xi_nonres_confirmed = np.array(xi_nonres_confirmed)

#To compute the K-S distances and their positions, as well as additional statistics:
delta_f = np.abs(len(P_obs)/(float(N_sim)/cos_factor) - len(P_confirmed)/float(N_Kep)) #absolute difference in the rates of observed planets per star

M_KS, M_KS_pos = KS_dist_mult(Mtot_obs[Mtot_obs > 0], M_confirmed)
P_KS, P_KS_pos = KS_dist(P_obs_flat, P_confirmed)
R_KS, R_KS_pos = KS_dist(Rm_obs_flat, R_confirmed)

R_res32_sim, R_res32_confirmed = np.float(sum((Rm_obs_flat >= res_ratios[0]) & (Rm_obs_flat <= res_ratios[0]*(1.+res_width))))/np.float(len(Rm_obs_flat)), np.float(sum((R_confirmed >= res_ratios[0]) & (R_confirmed <= res_ratios[0]*(1.+res_width))))/np.float(len(R_confirmed)) #fractions of planet pairs within 5% of 3:2 MMR, for simulated and Kepler data
R_res21_sim, R_res21_confirmed = np.float(sum((Rm_obs_flat >= res_ratios[1]) & (Rm_obs_flat <= res_ratios[1]*(1.+res_width))))/np.float(len(Rm_obs_flat)), np.float(sum((R_confirmed >= res_ratios[1]) & (R_confirmed <= res_ratios[1]*(1.+res_width))))/np.float(len(R_confirmed)) #fractions of planet pairs within 5% of 2:1 MMR, for simulated and Kepler data
R_res32_diff = np.abs(R_res32_sim - R_res32_confirmed) #difference in fractions of planet pairs close to 3:2 MMR between simulated and Kepler data
R_res21_diff = np.abs(R_res21_sim - R_res21_confirmed) #difference in fractions of planet pairs close to 2:1 MMR between simulated and Kepler data

tdur_KS, tdur_KS_pos = KS_dist(tdur_obs_flat, t_D_confirmed*60.)
D_KS, D_KS_pos = KS_dist(D_obs_flat, D_confirmed)
radii_KS, radii_KS_pos = KS_dist(radii_obs_flat, radii_confirmed)
D_ratio_KS, D_ratio_KS_pos = KS_dist(D_ratio_obs_flat, D_ratio_confirmed)
logxi_KS, logxi_KS_pos = KS_dist(np.log10(xi_obs_flat), np.log10(xi_confirmed))
logxi_res_KS, logxi_res_KS_pos = KS_dist(np.log10(xi_res_obs_flat), np.log10(xi_res_confirmed))
logxi_res32_KS, logxi_res32_KS_pos = KS_dist(np.log10(xi_res32_obs_flat), np.log10(xi_res32_confirmed))
logxi_res21_KS, logxi_res21_KS_pos = KS_dist(np.log10(xi_res21_obs_flat), np.log10(xi_res21_confirmed))
logxi_nonres_KS, logxi_nonres_KS_pos = KS_dist(np.log10(xi_nonres_obs_flat), np.log10(xi_nonres_confirmed))

distances = [delta_f, M_KS, P_KS, R_KS, tdur_KS, logxi_KS, D_KS, D_ratio_KS]

print 'Distances for (delta_f, Multiplicity, P, P ratio, t_dur, xi, depth, depth ratio):'
print 'Distances: ', [float(format(x, '.5f')) for x in distances]
print 'Total distance: ', sum(distances)



#To plot the 'observed' distributions with the actual observed Kepler distributions (observed total planet multiplicity, periods, period ratios):
'''
fig = plt.figure(figsize=(16,8))
plot = GridSpec(2,2,left=0.075,bottom=0.115,right=0.95,top=0.925,wspace=0.15,hspace=0.4)

ax = plt.subplot(plot[0,0])
plt.title('Observed planet multiplicity', fontsize=24)
x = Mtot_obs[Mtot_obs > 0]
max_M = np.max((np.max(Mtot_obs), np.max(M_confirmed)))
counts, bins = np.histogram(x, bins=max_M+1, range=(-0.5, max_M+0.5))
bins_mid = (bins[:-1] + bins[1:])/2.
plt.plot(bins_mid, counts/float(np.sum(counts)), 'o-', color='k', label='%s simulated observed systems' % len(x))
counts, bins = np.histogram(M_confirmed, bins=bins)
plt.plot(bins_mid, counts/float(np.sum(counts)), 'o--', color='k', alpha=0.2, label='%s Kepler systems' % len(M_confirmed))
ax.tick_params(axis='both', labelsize=20)
plt.xlim([1., max_M])
plt.xlabel(r'$M_{\rm tot}$', fontsize=20)
plt.ylabel('Fraction', fontsize=20)
plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[0,1])
plt.title('Observed periods', fontsize=24)
hist = plt.hist(P_obs_flat, bins=np.logspace(np.log10(np.min(P_obs_flat)), np.log10(np.max(P_obs_flat)), 21), histtype='step', weights=np.ones(len(P_obs_flat))/len(P_obs_flat), log=True, color='k', label=r'Simulated')
plt.hist(P_confirmed, bins=hist[1], histtype='stepfilled', weights=np.ones(len(P_confirmed))/len(P_confirmed), log=True, color='k', alpha=0.2, label=r'Kepler')
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([np.min(P_obs_flat)-1, 1.1*np.max(P_obs_flat)])
plt.ylim([10.**(-3.), 1.])
plt.xlabel(r'$P$ (days)', fontsize=20)
#plt.ylabel('Fraction', fontsize=20)
#plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[1,0])
plt.title('Observed period ratios ($P_{i+1}/P_i < 5$)', fontsize=24)
x = Rm_obs_flat[Rm_obs_flat < 5]
hist = plt.hist(x, bins=np.linspace(1,5,51), histtype='step', weights=np.ones(len(x))/len(x), color='k', label='Simulated')
x = R_confirmed[R_confirmed < 5]
plt.hist(x, bins=hist[1], histtype='stepfilled', weights=np.ones(len(x))/len(x), color='k', alpha=0.2, label='Kepler')
#plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([1,5])
plt.xlabel(r'$P_{i+1}/P_i$', fontsize=20)
plt.ylabel('Fraction', fontsize=20)
#plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[1,1])
plt.title('Observed period ratios ($P_{i+1}/P_i > 5$)', fontsize=24)
x = Rm_obs_flat[Rm_obs_flat > 5]
max_Rm = np.max((np.max(Rm_obs_flat), np.max(R_confirmed)))
hist = plt.hist(x, bins=np.linspace(5,max_Rm,51), histtype='step', weights=np.ones(len(x))/len(x), color='k', label='Simulated')
x = R_confirmed[R_confirmed > 5]
plt.hist(x, bins=hist[1], histtype='stepfilled', weights=np.ones(len(x))/len(x), color='k', alpha=0.2, label='Kepler')
#plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([5,max_Rm])
plt.xlabel(r'$P_{i+1}/P_i$', fontsize=20)
#plt.ylabel('Fraction', fontsize=20)
#plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

if savefigures == True:
    plt.savefig(savefigures_directory + model_name + '_observed_compare.pdf')
else:
    plt.show()
plt.close()
'''

#To plot the resulting 'observed' distributions with the actual observed Kepler distributions (observed transit durations, period-normalized transit duration ratios, transit depths, and transit depth ratios):
'''
fig = plt.figure(figsize=(16,8))
plot = GridSpec(2,2,left=0.075,bottom=0.115,right=0.95,top=0.925,wspace=0.15,hspace=0.4)

ax = plt.subplot(plot[0,0])
plt.title('Transit durations', fontsize=24)
x = tdur_obs_flat
hist = plt.hist(x, bins=50, histtype='step', weights=np.ones(len(x))/len(x), color='k', label='Simulated')
plt.hist(t_D_confirmed*60., bins=hist[1], histtype='stepfilled', weights=np.ones(len(t_D_confirmed))/len(t_D_confirmed), color='k', alpha=0.2, label=r'Kepler')
ax.tick_params(axis='both', labelsize=20)
plt.xlim([np.min(x), np.max(x)])
plt.xlabel(r'$t_{\rm dur}$ (mins)', fontsize=20)
plt.ylabel('Fraction', fontsize=20)

ax = plt.subplot(plot[1,0])
plt.title('Transit depths', fontsize=24)
x = D_obs_flat
hist = plt.hist(x, bins=np.logspace(np.log10(np.min(x)), np.log10(np.max(x)), 51), histtype='step', weights=np.ones(len(x))/len(x), color='k', label='Simulated')
plt.hist(D_confirmed, bins=hist[1], histtype='stepfilled', weights=np.ones(len(D_confirmed))/len(D_confirmed), color='k', alpha=0.2, label=r'Kepler')
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([np.min(x), np.max(x)])
plt.xlabel(r'$\delta$', fontsize=20)
plt.ylabel('Fraction', fontsize=20)
#plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[0,1])
plt.title('Period-normalized transit duration ratios', fontsize=24)
x = np.log10(xi_obs_flat)
hist = plt.hist(x, bins=np.linspace(-0.5,0.5,51), histtype='step', weights=np.ones(len(x))/len(x), color='k', label='All')
plt.hist(np.log10(xi_confirmed), bins=hist[1], histtype='stepfilled', weights=np.ones(len(xi_confirmed))/float(len(xi_confirmed)), color='k', alpha=0.2)
ax.tick_params(axis='both', labelsize=20)
plt.xlim([hist[1][0], hist[1][-1]])
plt.xlabel(r'$\log{\/xi}$', fontsize=20)
#plt.ylabel('Fraction', fontsize=20)
#plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[1,1])
plt.title('Transit depth ratios', fontsize=24)
x = D_ratio_obs_flat
hist = plt.hist(x, bins=np.logspace(np.log10(np.min(x)), np.log10(np.max(x)), 51), histtype='step', weights=np.ones(len(x))/len(x), color='k', label='Simulated')
plt.hist(D_ratio_confirmed, bins=hist[1], histtype='stepfilled', weights=np.ones(len(D_ratio_confirmed))/len(D_ratio_confirmed), color='k', alpha=0.2, label=r'Kepler')
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([hist[1][0], hist[1][-1]])
plt.xlabel(r'$\delta_{i+1}/\delta_i$', fontsize=20)
#plt.ylabel('Fraction', fontsize=20)
#plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

if savefigures == True:
    plt.savefig(savefigures_directory + model_name + '_observed_properties_compare.pdf')
else:
    plt.show()
plt.close()
'''

#To plot the resulting 'observed' distributions with the actual observed Kepler distributions (period-normalized transit duration ratios):
'''
fig = plt.figure(figsize=(16,8))
plot = GridSpec(2,2,left=0.075,bottom=0.115,right=0.95,top=0.925,wspace=0.15,hspace=0.4)

ax = plt.subplot(plot[0,0])
plt.title('Period-normalized transit duration ratios', fontsize=24)
x = np.log10(xi_obs_flat)
hist = plt.hist(x, bins=np.linspace(-0.5,0.5,51), histtype='step', weights=np.ones(len(x))/len(x), color='k', label='All')
plt.hist(np.log10(xi_confirmed), bins=hist[1], histtype='stepfilled', weights=np.ones(len(xi_confirmed))/float(len(xi_confirmed)), color='k', alpha=0.2)
plt.hist(np.log10(xi_res_obs_flat), bins=hist[1], histtype='step', weights=np.ones(len(xi_res_obs_flat))/float(len(x)), color='m', label='Near resonance')
plt.hist(np.log10(xi_res_confirmed), bins=hist[1], histtype='stepfilled', weights=np.ones(len(xi_res_confirmed))/float(len(xi_confirmed)), color='m', alpha=0.2)
plt.hist(np.log10(xi_nonres_obs_flat), bins=hist[1], histtype='step', weights=np.ones(len(xi_nonres_obs_flat))/float(len(xi_nonres_obs_flat)), color='g', label='Not near resonance')
plt.hist(np.log10(xi_nonres_confirmed), bins=hist[1], histtype='stepfilled', weights=np.ones(len(xi_nonres_confirmed))/float(len(xi_confirmed)), color='g', alpha=0.2)
ax.tick_params(axis='both', labelsize=20)
plt.xlim([hist[1][0], hist[1][-1]])
plt.xlabel(r'$\log{\/xi}$', fontsize=20)
plt.ylabel('Fraction', fontsize=20)
plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[1,0])
plt.title('Period-normalized transit duration ratios', fontsize=24)
x = np.log10(xi_obs_flat)
hist = plt.hist(x, bins=np.linspace(-0.5,0.5,51), histtype='step', weights=np.ones(len(x))/len(x), color='k', label='All')
plt.hist(np.log10(xi_res_obs_flat), bins=hist[1], histtype='step', weights=np.ones(len(xi_res_obs_flat))/float(len(xi_res_obs_flat)), color='m', label='Near resonance')
plt.hist(np.log10(xi_res_confirmed), bins=hist[1], histtype='stepfilled', weights=np.ones(len(xi_res_confirmed))/float(len(xi_res_confirmed)), color='m', alpha=0.2)
plt.hist(np.log10(xi_nonres_obs_flat), bins=hist[1], histtype='step', weights=np.ones(len(xi_nonres_obs_flat))/float(len(xi_nonres_obs_flat)), color='g', label='Not near resonance')
plt.hist(np.log10(xi_nonres_confirmed), bins=hist[1], histtype='stepfilled', weights=np.ones(len(xi_nonres_confirmed))/float(len(xi_nonres_confirmed)), color='g', alpha=0.2)
ax.tick_params(axis='both', labelsize=20)
plt.xlim([hist[1][0], hist[1][-1]])
plt.xlabel(r'$\log{\/xi}$', fontsize=20)
plt.ylabel('Fraction', fontsize=20)
plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[1,1])
plt.title('Period-normalized transit duration ratios', fontsize=24)
plt.hist(np.log10(xi_res32_obs_flat), bins=hist[1], histtype='step', weights=np.ones(len(xi_res32_obs_flat))/float(len(xi_res32_obs_flat)), color='r', label='Near 3:2 MMR')
plt.hist(np.log10(xi_res32_confirmed), bins=hist[1], histtype='stepfilled', weights=np.ones(len(xi_res32_confirmed))/float(len(xi_res32_confirmed)), color='r', alpha=0.2)
plt.hist(np.log10(xi_res21_obs_flat), bins=hist[1], histtype='step', weights=np.ones(len(xi_res21_obs_flat))/float(len(xi_res21_obs_flat)), color='b', label='Near 2:1 MMR')
plt.hist(np.log10(xi_res21_confirmed), bins=hist[1], histtype='stepfilled', weights=np.ones(len(xi_res21_confirmed))/float(len(xi_res21_confirmed)), color='b', alpha=0.2)
ax.tick_params(axis='both', labelsize=20)
plt.xlim([hist[1][0], hist[1][-1]])
plt.xlabel(r'$\log{\/xi}$', fontsize=20)
#plt.ylabel('Fraction', fontsize=20)
plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

if savefigures == True:
    plt.savefig(savefigures_directory + model_name + '_observed_logxi_compare.pdf')
else:
    plt.show()
plt.close()
'''


#To plot the CDFs (Kepler dataset and simulated observed systems) of the three observables:
'''
fig = plt.figure(figsize=(16,8))
plot = GridSpec(2,2,left=0.075,bottom=0.115,right=0.95,top=0.925,wspace=0.15,hspace=0.4)

ax = plt.subplot(plot[0,0])
plt.title('Observed planet multiplicity', fontsize=24)
x = Mtot_obs[Mtot_obs > 0]
counts_cumu = np.array([sum(x <= xi) for xi in range(1, np.max(x)+1)])
plt.plot(range(1, np.max(x)+1), counts_cumu/np.float(len(x)), drawstyle='steps-post', color='k', label='Simulated')
counts_cumu = np.array([sum(M_confirmed <= xi) for xi in range(1, np.max(M_confirmed)+1)])
plt.plot(range(1, np.max(M_confirmed)+1), counts_cumu/np.float(len(M_confirmed)), drawstyle='steps-post', color='k', ls='--', label='Kepler')
plt.annotate(s='', xy=(M_KS_pos + 0.5, sum(x <= M_KS_pos)/np.float(len(x))), xytext=(M_KS_pos + 0.5, sum(M_confirmed <= M_KS_pos)/np.float(len(M_confirmed))), arrowprops=dict(arrowstyle='<->'))
plt.text(np.max((np.max(x), np.max(M_confirmed))) + 0.8, 0.64, 'K-S distance = %s' % np.round(M_KS, 3), ha='right', fontsize=16)
ax.tick_params(axis='both', labelsize=20)
plt.xlim([1, np.max((np.max(x), np.max(M_confirmed)))+1])
plt.ylim([0.6,1])
plt.xlabel(r'$M_{\rm tot}$', fontsize=20)
plt.ylabel('Cumulative Fraction', fontsize=20)
plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[0,1])
plt.title('Observed periods $P_m$', fontsize=24)
plt.plot(np.sort(P_obs_flat), (np.arange(len(P_obs_flat))+1.)/np.float(len(P_obs_flat)), drawstyle='steps-post', color='k', label='Simulated')
plt.plot(np.sort(P_confirmed), (np.arange(len(P_confirmed))+1.)/np.float(len(P_confirmed)), drawstyle='steps-post', color='k', ls='--', label='Kepler')
plt.annotate(s='', xy=(P_KS_pos, sum(P_obs_flat <= P_KS_pos)/np.float(len(P_obs_flat))), xytext=(P_KS_pos, sum(P_confirmed <= P_KS_pos)/np.float(len(P_confirmed))), arrowprops=dict(arrowstyle='<->'))
plt.text(np.max(P_confirmed) - 100, 0.1, 'K-S distance = %s' % np.round(P_KS, 3), ha='right', fontsize=16)
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([0, np.max(P_confirmed)])
plt.xlabel(r'$P$ (days)', fontsize=20)
#plt.ylabel('Cumulative Fraction', fontsize=20)
#plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[1,0])
plt.title('Observed period ratios', fontsize=24)
plt.plot(np.sort(Rm_obs_flat), (np.arange(len(Rm_obs_flat))+1.)/np.float(len(Rm_obs_flat)), drawstyle='steps-post', color='k', label='Simulated')
plt.plot(np.sort(R_confirmed), (np.arange(len(R_confirmed))+1.)/np.float(len(R_confirmed)), drawstyle='steps-post', color='k', ls='--', label='Kepler')
plt.annotate(s='', xy=(R_KS_pos, sum(Rm_obs_flat <= R_KS_pos)/np.float(len(Rm_obs_flat))), xytext=(R_KS_pos, sum(R_confirmed <= R_KS_pos)/np.float(len(R_confirmed))), arrowprops=dict(arrowstyle='<->'))
plt.text(np.max((np.max(Rm_obs_flat), np.max(R_confirmed))) - 5, 0.1, 'K-S distance = %s' % np.round(R_KS, 3), ha='right', fontsize=16)
plt.text(np.max((np.max(Rm_obs_flat), np.max(R_confirmed))) - 5, 0.2, r'$\Delta{f_{3:2}} = |%s - %s| = %s$' % (np.round(R_res32_sim, 3), np.round(R_res32_confirmed, 3), np.round(R_res32_diff, 3)), ha='right', fontsize=16)
plt.text(np.max((np.max(Rm_obs_flat), np.max(R_confirmed))) - 5, 0.3, r'$\Delta{f_{2:1}} = |%s - %s| = %s$' % (np.round(R_res21_sim, 3), np.round(R_res21_confirmed, 3), np.round(R_res21_diff, 3)), ha='right', fontsize=16)
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([0, np.max((np.max(Rm_obs_flat), np.max(R_confirmed)))])
plt.xlabel(r'$P_{i+1}/P_i$', fontsize=20)
plt.ylabel('Cumulative Fraction', fontsize=20)
#plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

if savefigures == True:
    plt.savefig(savefigures_directory + model_name + '_KSdistances.pdf')
else:
    plt.show()
plt.close()
'''


#To plot the CDFs (Kepler dataset and simulated observed systems) of the transit durations, period-normalized transit duration ratios, transit depths, and transit depth ratios:
'''
fig = plt.figure(figsize=(16,8))
plot = GridSpec(2,2,left=0.075,bottom=0.115,right=0.95,top=0.925,wspace=0.15,hspace=0.4)

ax = plt.subplot(plot[0,0])
plt.title('Transit durations', fontsize=24)
plt.plot(np.sort(tdur_obs_flat), (np.arange(len(tdur_obs_flat))+1.)/np.float(len(tdur_obs_flat)), drawstyle='steps-post', color='k', label='Simulated')
plt.plot(np.sort(t_D_confirmed*60.), (np.arange(len(t_D_confirmed))+1.)/np.float(len(t_D_confirmed)), drawstyle='steps-post', color='k', ls='--', label='Kepler')
plt.annotate(s='', xy=(tdur_KS_pos, sum(tdur_obs_flat <= tdur_KS_pos)/np.float(len(tdur_obs_flat))), xytext=(tdur_KS_pos, sum(t_D_confirmed*60. <= tdur_KS_pos)/np.float(len(t_D_confirmed))), arrowprops=dict(arrowstyle='<->'))
plt.text(np.max((np.max(tdur_obs_flat), np.max(t_D_confirmed))) - 50, 0.1, 'K-S distance = %s' % np.round(tdur_KS, 3), ha='right', fontsize=16)
ax.tick_params(axis='both', labelsize=20)
plt.xlim([0, np.max(tdur_obs_flat)])
plt.xlabel(r'$t_{\rm dur}$ (mins)', fontsize=20)
plt.ylabel('Cumulative Fraction', fontsize=20)
plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[1,0])
plt.title('Transit depths', fontsize=24)
plt.plot(np.sort(D_obs_flat), (np.arange(len(D_obs_flat))+1.)/np.float(len(D_obs_flat)), drawstyle='steps-post', color='k', label='Simulated')
plt.plot(np.sort(D_confirmed), (np.arange(len(D_confirmed))+1.)/np.float(len(D_confirmed)), drawstyle='steps-post', color='k', ls='--', label='Kepler')
plt.annotate(s='', xy=(D_KS_pos, sum(D_obs_flat <= D_KS_pos)/np.float(len(D_obs_flat))), xytext=(D_KS_pos, sum(D_confirmed <= D_KS_pos)/np.float(len(D_confirmed))), arrowprops=dict(arrowstyle='<->'))
plt.text(np.max(D_obs_flat) - 0.005, 0.1, 'K-S distance = %s' % np.round(D_KS, 3), ha='right', fontsize=16)
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([0, np.max(D_obs_flat)])
plt.xlabel(r'$\delta$', fontsize=20)
plt.ylabel('Cumulative Fraction', fontsize=20)
#plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[0,1])
plt.title('Period-normalized transit duration ratios', fontsize=24)
plt.plot(np.sort(np.log10(xi_obs_flat)), (np.arange(len(np.log10(xi_obs_flat)))+1.)/np.float(len(np.log10(xi_obs_flat))), drawstyle='steps-post', color='k', label='All')
plt.plot(np.sort(np.log10(xi_confirmed)), (np.arange(len(np.log10(xi_confirmed)))+1.)/np.float(len(np.log10(xi_confirmed))), drawstyle='steps-post', color='k', ls='--')
plt.annotate(s='', xy=(logxi_KS_pos, sum(np.log10(xi_obs_flat) <= logxi_KS_pos)/np.float(len(xi_obs_flat))), xytext=(logxi_KS_pos, sum(np.log10(xi_confirmed) <= logxi_KS_pos)/np.float(len(xi_confirmed))), arrowprops=dict(arrowstyle='<->'))
plt.text(0.45, 0.1, 'K-S distance = %s' % np.round(logxi_KS, 3), ha='right', fontsize=16)
ax.tick_params(axis='both', labelsize=20)
plt.xlim([-0.5, 0.5])
plt.xlabel(r'$\log{\/xi}$', fontsize=20)
#plt.ylabel('Cumulative Fraction', fontsize=20)
#plt.legend(loc='upper left', bbox_to_anchor=(0.05,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[1,1])
plt.title('Transit depth ratios', fontsize=24)
plt.plot(np.sort(D_ratio_obs_flat), (np.arange(len(D_ratio_obs_flat))+1.)/np.float(len(D_ratio_obs_flat)), drawstyle='steps-post', color='k', label='All')
plt.plot(np.sort(D_ratio_confirmed), (np.arange(len(D_ratio_confirmed))+1.)/np.float(len(D_ratio_confirmed)), drawstyle='steps-post', color='k', ls='--')
plt.annotate(s='', xy=(D_ratio_KS_pos, sum(D_ratio_obs_flat <= D_ratio_KS_pos)/np.float(len(D_ratio_obs_flat))), xytext=(D_ratio_KS_pos, sum(D_ratio_confirmed <= D_ratio_KS_pos)/np.float(len(D_ratio_confirmed))), arrowprops=dict(arrowstyle='<->'))
plt.text(0.8*np.max((np.max(D_ratio_obs_flat), np.max(D_ratio_confirmed))), 0.1, 'K-S distance = %s' % np.round(D_ratio_KS, 3), ha='right', fontsize=16)
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([0, np.max((np.max(D_ratio_obs_flat), np.max(D_ratio_confirmed)))])
plt.xlabel(r'$\delta_{i+1}/\delta_i$', fontsize=20)
#plt.ylabel('Cumulative Fraction', fontsize=20)
#plt.legend(loc='upper left', bbox_to_anchor=(0.05,0.95), ncol=1, fontsize=16) #show the legend

if savefigures == True:
    plt.savefig(savefigures_directory + model_name + '_KSdistances_properties.pdf')
else:
    plt.show()
plt.close()
'''

#To plot the CDFs (Kepler dataset and simulated observed systems) of the period-normalized transit duration ratios:
'''
fig = plt.figure(figsize=(16,8))
plot = GridSpec(2,2,left=0.075,bottom=0.115,right=0.95,top=0.925,wspace=0.15,hspace=0.4)

ax = plt.subplot(plot[0,0])
plt.title('Period-normalized transit duration ratios', fontsize=24)
plt.plot(np.sort(np.log10(xi_obs_flat)), (np.arange(len(np.log10(xi_obs_flat)))+1.)/np.float(len(np.log10(xi_obs_flat))), drawstyle='steps-post', color='k', label='All')
plt.plot(np.sort(np.log10(xi_confirmed)), (np.arange(len(np.log10(xi_confirmed)))+1.)/np.float(len(np.log10(xi_confirmed))), drawstyle='steps-post', color='k', ls='--')
plt.text(0.45, 0.1, 'K-S distance = %s' % np.round(logxi_KS, 3), ha='right', fontsize=16)
ax.tick_params(axis='both', labelsize=20)
plt.xlim([-0.5, 0.5])
plt.xlabel(r'$\log{\/xi}$', fontsize=20)
#plt.ylabel('Cumulative Fraction', fontsize=20)
plt.legend(loc='upper left', bbox_to_anchor=(0.05,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[1,0])
plt.title('Period-normalized transit duration ratios', fontsize=24)
plt.plot(np.sort(np.log10(xi_res_obs_flat)), (np.arange(len(np.log10(xi_res_obs_flat)))+1.)/np.float(len(np.log10(xi_res_obs_flat))), drawstyle='steps-post', color='m', label='Near resonance')
plt.plot(np.sort(np.log10(xi_res_confirmed)), (np.arange(len(np.log10(xi_res_confirmed)))+1.)/np.float(len(np.log10(xi_res_confirmed))), drawstyle='steps-post', color='m', ls='--')
plt.plot(np.sort(np.log10(xi_nonres_obs_flat)), (np.arange(len(np.log10(xi_nonres_obs_flat)))+1.)/np.float(len(np.log10(xi_nonres_obs_flat))), drawstyle='steps-post', color='g', label='Not near resonance')
plt.plot(np.sort(np.log10(xi_nonres_confirmed)), (np.arange(len(np.log10(xi_nonres_confirmed)))+1.)/np.float(len(np.log10(xi_nonres_confirmed))), drawstyle='steps-post', color='g', ls='--')
plt.text(0.45, 0.2, 'K-S distance = %s' % np.round(logxi_res_KS, 3), ha='right', fontsize=16)
plt.text(0.45, 0.1, 'K-S distance = %s' % np.round(logxi_nonres_KS, 3), ha='right', fontsize=16)
ax.tick_params(axis='both', labelsize=20)
plt.xlim([-0.5, 0.5])
plt.xlabel(r'$\log{\/xi}$', fontsize=20)
#plt.ylabel('Cumulative Fraction', fontsize=20)
plt.legend(loc='upper left', bbox_to_anchor=(0.05,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[1,1])
plt.title('Period-normalized transit duration ratios', fontsize=24)
plt.plot(np.sort(np.log10(xi_res32_obs_flat)), (np.arange(len(np.log10(xi_res32_obs_flat)))+1.)/np.float(len(np.log10(xi_res32_obs_flat))), drawstyle='steps-post', color='r', label='Near 3:2 MMR')
plt.plot(np.sort(np.log10(xi_res32_confirmed)), (np.arange(len(np.log10(xi_res32_confirmed)))+1.)/np.float(len(np.log10(xi_res32_confirmed))), drawstyle='steps-post', color='r', ls='--')
plt.plot(np.sort(np.log10(xi_res21_obs_flat)), (np.arange(len(np.log10(xi_res21_obs_flat)))+1.)/np.float(len(np.log10(xi_res21_obs_flat))), drawstyle='steps-post', color='b', label='Near 2:1 MMR')
plt.plot(np.sort(np.log10(xi_res21_confirmed)), (np.arange(len(np.log10(xi_res21_confirmed)))+1.)/np.float(len(np.log10(xi_res21_confirmed))), drawstyle='steps-post', color='b', ls='--')
plt.text(0.45, 0.2, 'K-S distance = %s' % np.round(logxi_res32_KS, 3), ha='right', fontsize=16)
plt.text(0.45, 0.1, 'K-S distance = %s' % np.round(logxi_res21_KS, 3), ha='right', fontsize=16)
ax.tick_params(axis='both', labelsize=20)
plt.xlim([-0.5, 0.5])
plt.xlabel(r'$\log{\/xi}$', fontsize=20)
#plt.ylabel('Cumulative Fraction', fontsize=20)
plt.legend(loc='upper left', bbox_to_anchor=(0.05,0.95), ncol=1, fontsize=16) #show the legend

if savefigures == True:
    plt.savefig(savefigures_directory + model_name + '_KSdistances_logxi.pdf')
else:
    plt.show()
plt.close()
'''





#'''
##### To re-make all those comparison and CDF plots as single panels for 2nd year paper or presentation:

subdirectory = 'Talk_Figures/' #'Paper_Figures/'; 'Talk_Figures/'

if subdirectory == 'Paper_Figures/':
    linewidth = 1
elif subdirectory == 'Talk_Figures/':
    linewidth = 3

#To make a 'plot' listing the model parameters:
fig = plt.figure(figsize=(8,4))
plot = GridSpec(1,1,left=0.15,bottom=0.2,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
nrows = 8
for i in range(len(param_keys_all)): #range(len(param_keys_all))
    plt.figtext(x=0.05+0.3*int(i/float(nrows)), y=0.875-0.1*(i%nrows), s=r'%s = %s' % (param_keys_all[i][1], param_vals_all[i]), fontsize=16)
if savefigures == True:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_sim_params.pdf')
    plt.close()

#To plot the multiplicities:
fig = plt.figure(figsize=(8,4))
plot = GridSpec(1,1,left=0.15,bottom=0.2,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
x = Mtot_obs[Mtot_obs > 0]
max_M = np.max((np.max(Mtot_obs), np.max(M_confirmed)))
counts, bins = np.histogram(x, bins=max_M+1, range=(-0.5, max_M+0.5))
bins_mid = (bins[:-1] + bins[1:])/2.
plt.plot(bins_mid, counts/float(np.sum(counts)), 'o-', color='k', linewidth=linewidth, label='%s simulated observed systems' % len(x))
counts, bins = np.histogram(M_confirmed, bins=bins)
plt.plot(bins_mid, counts/float(np.sum(counts)), 'o--', color='k', label='%s Kepler systems' % len(M_confirmed))
plt.gca().set_yscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([1., max_M])
plt.xlabel(r'Number of observed planets', fontsize=20)
plt.ylabel('Fraction', fontsize=20)
plt.legend(loc='upper right', bbox_to_anchor=(0.99,0.99), ncol=1, frameon=False, fontsize=16) #show the legend
if savefigures == True:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_multiplicities_compare.pdf')
    plt.close()

#To plot the periods:
fig = plt.figure(figsize=(8,4))
plot = GridSpec(1,1,left=0.15,bottom=0.2,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
hist = plt.hist(P_obs_flat, bins=np.logspace(np.log10(np.min(P_obs_flat)), np.log10(np.max(P_obs_flat)), 101), histtype='step', weights=np.ones(len(P_obs_flat))/len(P_obs_flat), log=True, color='k', linewidth=linewidth, label=r'Simulated')
plt.hist(P_confirmed, bins=hist[1], histtype='stepfilled', weights=np.ones(len(P_confirmed))/len(P_confirmed), log=True, color='k', alpha=0.2, label=r'Kepler')
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
ax.set_xticks([3,10,30,100,300])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.xlim([P_min, 1.1*P_max])
plt.ylim([10.**(-3.), 0.1])
plt.xlabel(r'$P$ (days)', fontsize=20)
plt.ylabel('Fraction', fontsize=20)
plt.legend(loc='upper right', bbox_to_anchor=(0.99,0.99), ncol=1, frameon=False, fontsize=16) #show the legend
if savefigures == True:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_periods_compare.pdf')
    plt.close()

#To plot the period ratios (all):
fig = plt.figure(figsize=(8,4))
plot = GridSpec(1,1,left=0.15,bottom=0.2,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
R_max_cut = 30. #upper cut-off for plotting period ratios; np.max(Rm_obs_flat)
x = Rm_obs_flat[Rm_obs_flat < R_max_cut]
hist = plt.hist(x, bins=np.logspace(np.log10(min(np.min(x),np.min(R_confirmed))), np.log10(R_max_cut), 101), histtype='step', weights=np.ones(len(x))/len(x), color='k', linewidth=linewidth, label='Simulated')
x = R_confirmed[R_confirmed < R_max_cut]
plt.hist(x, bins=hist[1], histtype='stepfilled', weights=np.ones(len(x))/len(x), color='k', alpha=0.2, label='Kepler')
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
ax.set_xticks([1,2,3,4,5,10,20])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.xlim([1, R_max_cut])
plt.xlabel(r'$P_{i+1}/P_i$', fontsize=20)
plt.ylabel('Fraction', fontsize=20)
if subdirectory == 'Talk_Figures/':
    plt.legend(loc='upper right', bbox_to_anchor=(0.99,0.99), ncol=1, frameon=False, fontsize=16) #show the legend
if savefigures == True:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_periodratios_compare.pdf')
    plt.close()

#To plot the period ratios (< 5):
fig = plt.figure(figsize=(8,4))
plot = GridSpec(1,1,left=0.15,bottom=0.2,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
x = Rm_obs_flat[Rm_obs_flat < 5]
hist = plt.hist(x, bins=np.linspace(1,5,101), histtype='step', weights=np.ones(len(x))/len(x), color='k', linewidth=linewidth, label='Simulated')
x = R_confirmed[R_confirmed < 5]
plt.hist(x, bins=hist[1], histtype='stepfilled', weights=np.ones(len(x))/len(x), color='k', alpha=0.2, label='Kepler')
ax.tick_params(axis='both', labelsize=20)
plt.xlim([1,5])
plt.xlabel(r'$P_{i+1}/P_i$', fontsize=20)
plt.ylabel('Fraction', fontsize=20)
if subdirectory == 'Talk_Figures/':
    plt.legend(loc='upper right', bbox_to_anchor=(0.99,0.99), ncol=1, frameon=False, fontsize=16) #show the legend
if savefigures == True:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_periodratios_less5_compare.pdf')
    plt.close()

#To plot the period ratios (> 5):
fig = plt.figure(figsize=(8,4))
plot = GridSpec(1,1,left=0.15,bottom=0.2,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
x = Rm_obs_flat[Rm_obs_flat > 5]
max_Rm = np.max((np.max(Rm_obs_flat), np.max(R_confirmed)))
hist = plt.hist(x, bins=np.linspace(5,max_Rm,101), histtype='step', weights=np.ones(len(x))/len(x), color='k', linewidth=linewidth, label='Simulated')
x = R_confirmed[R_confirmed > 5]
plt.hist(x, bins=hist[1], histtype='stepfilled', weights=np.ones(len(x))/len(x), color='k', alpha=0.2, label='Kepler')
ax.tick_params(axis='both', labelsize=20)
plt.xlim([5,np.max(x)])
plt.xlabel(r'$P_{i+1}/P_i$', fontsize=20)
plt.ylabel('Fraction', fontsize=20)
if subdirectory == 'Talk_Figures/':
    plt.legend(loc='upper right', bbox_to_anchor=(0.99,0.99), ncol=1, frameon=False, fontsize=16) #show the legend
if savefigures == True:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_periodratios_greater5_compare.pdf')
    plt.close()

#To plot the transit durations:
fig = plt.figure(figsize=(8,4))
plot = GridSpec(1,1,left=0.15,bottom=0.2,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
x = tdur_obs_flat
hist = plt.hist(x, bins=np.linspace(min(np.min(tdur_obs_flat),np.min(t_D_confirmed*60.)), max(np.max(tdur_obs_flat),np.max(t_D_confirmed*60.)), 101), histtype='step', weights=np.ones(len(x))/len(x), color='k', linewidth=linewidth, label='Simulated')
plt.hist(t_D_confirmed*60., bins=hist[1], histtype='stepfilled', weights=np.ones(len(t_D_confirmed))/len(t_D_confirmed), color='k', alpha=0.2, label=r'Kepler')
ax.tick_params(axis='both', labelsize=20)
plt.xlim([hist[1][0], hist[1][-1]])
plt.ylim([0., 1.1*np.max(hist[0])])
plt.xlabel(r'$t_{\rm dur}$ (mins)', fontsize=20)
plt.ylabel('Fraction', fontsize=20)
if subdirectory == 'Talk_Figures/':
    plt.legend(loc='upper right', bbox_to_anchor=(0.99,0.99), ncol=1, frameon=False, fontsize=16) #show the legend
if savefigures == True:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_durations_compare.pdf')
    plt.close()

#To plot the transit depths:
fig = plt.figure(figsize=(8,4))
plot = GridSpec(1,1,left=0.15,bottom=0.2,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
x = D_obs_flat
hist = plt.hist(x, bins=np.logspace(np.log10(min(np.min(x),np.min(D_confirmed))), np.log10(max(np.max(x),np.max(D_confirmed))), 101), histtype='step', weights=np.ones(len(x))/len(x), color='k', linewidth=linewidth, label='Simulated')
plt.hist(D_confirmed, bins=hist[1], histtype='stepfilled', weights=np.ones(len(D_confirmed))/len(D_confirmed), color='k', alpha=0.2, label=r'Kepler')
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([hist[1][0], hist[1][-1]])
plt.xlabel(r'$\delta$', fontsize=20)
plt.ylabel('Fraction', fontsize=20)
if subdirectory == 'Talk_Figures/':
    plt.legend(loc='upper right', bbox_to_anchor=(0.99,0.99), ncol=1, frameon=False, fontsize=16) #show the legend
if savefigures == True:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_depths_compare.pdf')
    plt.close()

#To plot the transit depth ratios:
fig = plt.figure(figsize=(8,4))
plot = GridSpec(1,1,left=0.15,bottom=0.2,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
x = D_ratio_obs_flat
hist = plt.hist(x, bins=np.logspace(np.log10(min(np.min(x),np.min(D_ratio_confirmed))), np.log10(max(np.max(x),np.max(D_ratio_confirmed))), 101), histtype='step', weights=np.ones(len(x))/len(x), color='k', linewidth=linewidth, label='Simulated')
plt.hist(D_ratio_confirmed, bins=hist[1], histtype='stepfilled', weights=np.ones(len(D_ratio_confirmed))/len(D_ratio_confirmed), color='k', alpha=0.2, label=r'Kepler')
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([hist[1][0], hist[1][-1]])
plt.xlabel(r'$\delta_{i+1}/\delta_i$', fontsize=20)
plt.ylabel('Fraction', fontsize=20)
if subdirectory == 'Talk_Figures/':
    plt.legend(loc='upper right', bbox_to_anchor=(0.99,0.99), ncol=1, frameon=False, fontsize=16) #show the legend
if savefigures == True:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_depthratios_compare.pdf')
    plt.close()

#To plot log(xi):
fig = plt.figure(figsize=(8,4))
plot = GridSpec(1,1,left=0.15,bottom=0.2,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
x = np.log10(xi_obs_flat)
hist = plt.hist(x, bins=np.linspace(-0.5, 0.5, 101), histtype='step', weights=np.ones(len(x))/len(x), color='k', linewidth=linewidth, label='All')
plt.hist(np.log10(xi_confirmed), bins=hist[1], histtype='stepfilled', weights=np.ones(len(xi_confirmed))/float(len(xi_confirmed)), color='k', alpha=0.2)
ax.tick_params(axis='both', labelsize=20)
plt.xlim([hist[1][0], hist[1][-1]])
plt.xlabel(r'$\log{\xi}$', fontsize=20)
plt.ylabel('Fraction', fontsize=20)
if subdirectory == 'Talk_Figures/':
    plt.legend(loc='upper right', bbox_to_anchor=(0.99,0.99), ncol=1, frameon=False, fontsize=16) #show the legend
if savefigures == True:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_logxi_all_compare.pdf')
    plt.close()

#To plot log(xi) by res/non-res:
fig = plt.figure(figsize=(8,4))
plot = GridSpec(1,1,left=0.15,bottom=0.2,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
x = np.log10(xi_obs_flat)
hist = plt.hist(x, bins=np.linspace(-0.5,0.5,101), histtype='step', weights=np.ones(len(x))/len(x), color='k', label='All')
plt.hist(np.log10(xi_res_obs_flat), bins=hist[1], histtype='step', weights=np.ones(len(xi_res_obs_flat))/float(len(xi_res_obs_flat)), color='m', linewidth=linewidth, label='Near resonance')
plt.hist(np.log10(xi_res_confirmed), bins=hist[1], histtype='stepfilled', weights=np.ones(len(xi_res_confirmed))/float(len(xi_res_confirmed)), color='m', alpha=0.2)
plt.hist(np.log10(xi_nonres_obs_flat), bins=hist[1], histtype='step', weights=np.ones(len(xi_nonres_obs_flat))/float(len(xi_nonres_obs_flat)), color='g', linewidth=linewidth, label='Not near resonance')
plt.hist(np.log10(xi_nonres_confirmed), bins=hist[1], histtype='stepfilled', weights=np.ones(len(xi_nonres_confirmed))/float(len(xi_nonres_confirmed)), color='g', alpha=0.2)
ax.tick_params(axis='both', labelsize=20)
plt.xlim([hist[1][0], hist[1][-1]])
plt.xlabel(r'$\log{\xi}$', fontsize=20)
plt.ylabel('Fraction', fontsize=20)
plt.legend(loc='upper right', bbox_to_anchor=(0.99,0.99), ncol=1, frameon=False, fontsize=16) #show the legend
if savefigures == True:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_logxi_compare.pdf')
    plt.close()

#To plot log(xi) within res:
fig = plt.figure(figsize=(8,4))
plot = GridSpec(1,1,left=0.15,bottom=0.2,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
plt.hist(np.log10(xi_res32_obs_flat), bins=hist[1], histtype='step', weights=np.ones(len(xi_res32_obs_flat))/float(len(xi_res32_obs_flat)), color='r', linewidth=linewidth, label='Near 3:2 MMR')
plt.hist(np.log10(xi_res32_confirmed), bins=hist[1], histtype='stepfilled', weights=np.ones(len(xi_res32_confirmed))/float(len(xi_res32_confirmed)), color='r', alpha=0.2)
plt.hist(np.log10(xi_res21_obs_flat), bins=hist[1], histtype='step', weights=np.ones(len(xi_res21_obs_flat))/float(len(xi_res21_obs_flat)), color='b', linewidth=linewidth, label='Near 2:1 MMR')
plt.hist(np.log10(xi_res21_confirmed), bins=hist[1], histtype='stepfilled', weights=np.ones(len(xi_res21_confirmed))/float(len(xi_res21_confirmed)), color='b', alpha=0.2)
ax.tick_params(axis='both', labelsize=20)
plt.xlim([hist[1][0], hist[1][-1]])
plt.xlabel(r'$\log{\xi}$', fontsize=20)
plt.ylabel('Fraction', fontsize=20)
plt.legend(loc='upper right', bbox_to_anchor=(0.99,0.99), ncol=1, frameon=False, fontsize=16) #show the legend
if savefigures == True:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_logxi_res_compare.pdf')
    plt.close()

plt.show()
plt.close()



#To plot the multiplicity CDFs:
fig = plt.figure(figsize=(8,3))
plot = GridSpec(1,1,left=0.15,bottom=0.3,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
x = Mtot_obs[Mtot_obs > 0]
counts_cumu = np.array([sum(x <= xi) for xi in range(1, np.max(x)+1)])
plt.plot(range(1, np.max(x)+1), counts_cumu/np.float(len(x)), drawstyle='steps-post', color='k', linewidth=linewidth, label='Simulated')
counts_cumu = np.array([sum(M_confirmed <= xi) for xi in range(1, np.max(M_confirmed)+1)])
plt.plot(range(1, np.max(M_confirmed)+1), counts_cumu/np.float(len(M_confirmed)), drawstyle='steps-post', color='k', ls='--', label='Kepler')
plt.annotate(s='', xy=(M_KS_pos + 0.5, sum(x <= M_KS_pos)/np.float(len(x))), xytext=(M_KS_pos + 0.5, sum(M_confirmed <= M_KS_pos)/np.float(len(M_confirmed))), arrowprops=dict(arrowstyle='<->'))
plt.figtext(0.925, 0.35, r'$\mathcal{D} = %s$' % np.round(M_KS, 3), ha='right', fontsize=16)
ax.tick_params(axis='both', labelsize=20)
plt.xlim([1, np.max((np.max(x), np.max(M_confirmed)))+1])
plt.ylim([0.6,1])
plt.xlabel(r'Number of planets', fontsize=20)
plt.ylabel('CDF', fontsize=20)
plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, frameon=False, fontsize=16) #show the legend
if savefigures == True:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_multiplicities_CDFs.pdf')
    plt.close()

#To plot the period CDFs:
fig = plt.figure(figsize=(8,3))
plot = GridSpec(1,1,left=0.15,bottom=0.3,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
plt.plot(np.sort(P_obs_flat), (np.arange(len(P_obs_flat))+1.)/np.float(len(P_obs_flat)), drawstyle='steps-post', color='k', linewidth=linewidth, label='Simulated')
plt.plot(np.sort(P_confirmed), (np.arange(len(P_confirmed))+1.)/np.float(len(P_confirmed)), drawstyle='steps-post', color='k', ls='--', label='Kepler')
plt.annotate(s='', xy=(P_KS_pos, sum(P_obs_flat <= P_KS_pos)/np.float(len(P_obs_flat))), xytext=(P_KS_pos, sum(P_confirmed <= P_KS_pos)/np.float(len(P_confirmed))), arrowprops=dict(arrowstyle='<->'))
plt.figtext(0.925, 0.35, r'$\mathcal{D} = %s$' % np.round(P_KS, 3), ha='right', fontsize=16)
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
ax.set_xticks([3,10,30,100,300])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.xlim([0, np.max(P_confirmed)])
plt.xlabel(r'$P$ (days)', fontsize=20)
plt.ylabel('CDF', fontsize=20)
if savefigures == True:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_periods_CDFs.pdf')
    plt.close()

#To plot the period ratio CDFs:
fig = plt.figure(figsize=(8,3))
plot = GridSpec(1,1,left=0.15,bottom=0.3,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
plt.plot(np.sort(Rm_obs_flat), (np.arange(len(Rm_obs_flat))+1.)/np.float(len(Rm_obs_flat)), drawstyle='steps-post', color='k', linewidth=linewidth, label='Simulated')
plt.plot(np.sort(R_confirmed), (np.arange(len(R_confirmed))+1.)/np.float(len(R_confirmed)), drawstyle='steps-post', color='k', ls='--', label='Kepler')
plt.annotate(s='', xy=(R_KS_pos, sum(Rm_obs_flat <= R_KS_pos)/np.float(len(Rm_obs_flat))), xytext=(R_KS_pos, sum(R_confirmed <= R_KS_pos)/np.float(len(R_confirmed))), arrowprops=dict(arrowstyle='<->'))
plt.figtext(0.925, 0.35, r'$\mathcal{D} = %s$' % np.round(R_KS, 3), ha='right', fontsize=16)
#plt.figtext(0.925, 0.4, r'$\Delta{f_{3:2}} = |%s - %s| = %s$' % (np.round(R_res32_sim, 3), np.round(R_res32_confirmed, 3), np.round(R_res32_diff, 3)), ha='right', fontsize=16)
#plt.figtext(0.925, 0.45, r'$\Delta{f_{2:1}} = |%s - %s| = %s$' % (np.round(R_res21_sim, 3), np.round(R_res21_confirmed, 3), np.round(R_res21_diff, 3)), ha='right', fontsize=16)
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
ax.set_xticks([1,2,3,4,5,10,20,40])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.xlim([0, np.max((np.max(Rm_obs_flat), np.max(R_confirmed)))])
plt.xlabel(r'$P_{i+1}/P_i$', fontsize=20)
plt.ylabel('CDF', fontsize=20)
if savefigures == True:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_periodratios_CDFs.pdf')
    plt.close()

#To plot the transit duration CDFs:
fig = plt.figure(figsize=(8,3))
plot = GridSpec(1,1,left=0.15,bottom=0.3,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
plt.plot(np.sort(tdur_obs_flat), (np.arange(len(tdur_obs_flat))+1.)/np.float(len(tdur_obs_flat)), drawstyle='steps-post', color='k', linewidth=linewidth, label='Simulated')
plt.plot(np.sort(t_D_confirmed*60.), (np.arange(len(t_D_confirmed))+1.)/np.float(len(t_D_confirmed)), drawstyle='steps-post', color='k', ls='--', label='Kepler')
plt.annotate(s='', xy=(tdur_KS_pos, sum(tdur_obs_flat <= tdur_KS_pos)/np.float(len(tdur_obs_flat))), xytext=(tdur_KS_pos, sum(t_D_confirmed*60. <= tdur_KS_pos)/np.float(len(t_D_confirmed))), arrowprops=dict(arrowstyle='<->'))
plt.figtext(0.925, 0.35, r'$\mathcal{D} = %s$' % np.round(tdur_KS, 3), ha='right', fontsize=16)
ax.tick_params(axis='both', labelsize=20)
plt.xlim([0, np.max(tdur_obs_flat)])
plt.xlabel(r'$t_{\rm dur}$ (mins)', fontsize=20)
plt.ylabel('CDF', fontsize=20)
if savefigures == True:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_durations_CDFs.pdf')
    plt.close()

#To plot the transit depth CDFs:
fig = plt.figure(figsize=(8,3))
plot = GridSpec(1,1,left=0.15,bottom=0.3,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
plt.plot(np.sort(D_obs_flat), (np.arange(len(D_obs_flat))+1.)/np.float(len(D_obs_flat)), drawstyle='steps-post', color='k', linewidth=linewidth, label='Simulated')
plt.plot(np.sort(D_confirmed), (np.arange(len(D_confirmed))+1.)/np.float(len(D_confirmed)), drawstyle='steps-post', color='k', ls='--', label='Kepler')
plt.annotate(s='', xy=(D_KS_pos, sum(D_obs_flat <= D_KS_pos)/np.float(len(D_obs_flat))), xytext=(D_KS_pos, sum(D_confirmed <= D_KS_pos)/np.float(len(D_confirmed))), arrowprops=dict(arrowstyle='<->'))
plt.figtext(0.925, 0.35, r'$\mathcal{D} = %s$' % np.round(D_KS, 3), ha='right', fontsize=16)
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([0, np.max(D_obs_flat)])
plt.xlabel(r'$\delta$', fontsize=20)
plt.ylabel('CDF', fontsize=20)
if savefigures == True:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_depths_CDFs.pdf')
    plt.close()

#To plot the transit depth ratio CDFs:
fig = plt.figure(figsize=(8,3))
plot = GridSpec(1,1,left=0.15,bottom=0.3,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
plt.plot(np.sort(D_ratio_obs_flat), (np.arange(len(D_ratio_obs_flat))+1.)/np.float(len(D_ratio_obs_flat)), drawstyle='steps-post', color='k', linewidth=linewidth, label='Simulated')
plt.plot(np.sort(D_ratio_confirmed), (np.arange(len(D_ratio_confirmed))+1.)/np.float(len(D_ratio_confirmed)), drawstyle='steps-post', color='k', ls='--', label='Kepler')
plt.annotate(s='', xy=(D_ratio_KS_pos, sum(D_ratio_obs_flat <= D_ratio_KS_pos)/np.float(len(D_ratio_obs_flat))), xytext=(D_ratio_KS_pos, sum(D_ratio_confirmed <= D_ratio_KS_pos)/np.float(len(D_ratio_confirmed))), arrowprops=dict(arrowstyle='<->'))
plt.figtext(0.925, 0.35, r'$\mathcal{D} = %s$' % np.round(D_ratio_KS, 3), ha='right', fontsize=16)
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([0, np.max((np.max(D_ratio_obs_flat), np.max(D_ratio_confirmed)))])
plt.xlabel(r'$\delta_{i+1}/\delta_i$', fontsize=20)
plt.ylabel('CDF', fontsize=20)
if savefigures == True:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_depthratios_CDFs.pdf')
    plt.close()

#To plot the log(xi) CDFs:
fig = plt.figure(figsize=(8,3))
plot = GridSpec(1,1,left=0.15,bottom=0.3,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
plt.plot(np.sort(np.log10(xi_obs_flat)), (np.arange(len(np.log10(xi_obs_flat)))+1.)/np.float(len(np.log10(xi_obs_flat))), drawstyle='steps-post', color='k', linewidth=linewidth, label='All')
plt.plot(np.sort(np.log10(xi_confirmed)), (np.arange(len(np.log10(xi_confirmed)))+1.)/np.float(len(np.log10(xi_confirmed))), drawstyle='steps-post', color='k', ls='--')
plt.annotate(s='', xy=(P_KS_pos, sum(P_obs_flat <= P_KS_pos)/np.float(len(P_obs_flat))), xytext=(P_KS_pos, sum(P_confirmed <= P_KS_pos)/np.float(len(P_confirmed))), arrowprops=dict(arrowstyle='<->'))
plt.figtext(0.925, 0.35, r'$\mathcal{D} = %s$' % np.round(logxi_KS, 3), ha='right', fontsize=16)
ax.tick_params(axis='both', labelsize=20)
plt.xlim([-0.5, 0.5])
plt.xlabel(r'$\log{\xi}$', fontsize=20)
plt.ylabel('CDF', fontsize=20)
if savefigures == True:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_logxi_all_CDFs.pdf')
    plt.close()

#To plot the log(xi) CDFs by res/non-res:
fig = plt.figure(figsize=(8,4))
plot = GridSpec(1,1,left=0.15,bottom=0.2,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
plt.plot(np.sort(np.log10(xi_obs_flat)), (np.arange(len(np.log10(xi_obs_flat)))+1.)/np.float(len(np.log10(xi_obs_flat))), drawstyle='steps-post', color='k', linewidth=linewidth, label='All')
plt.plot(np.sort(np.log10(xi_confirmed)), (np.arange(len(np.log10(xi_confirmed)))+1.)/np.float(len(np.log10(xi_confirmed))), drawstyle='steps-post', color='k', ls='--')
plt.plot(np.sort(np.log10(xi_res_obs_flat)), (np.arange(len(np.log10(xi_res_obs_flat)))+1.)/np.float(len(np.log10(xi_res_obs_flat))), drawstyle='steps-post', color='m', linewidth=linewidth, label='Near resonance')
plt.plot(np.sort(np.log10(xi_res_confirmed)), (np.arange(len(np.log10(xi_res_confirmed)))+1.)/np.float(len(np.log10(xi_res_confirmed))), drawstyle='steps-post', color='m', ls='--')
plt.plot(np.sort(np.log10(xi_nonres_obs_flat)), (np.arange(len(np.log10(xi_nonres_obs_flat)))+1.)/np.float(len(np.log10(xi_nonres_obs_flat))), drawstyle='steps-post', color='g', linewidth=linewidth, label='Not near resonance')
plt.plot(np.sort(np.log10(xi_nonres_confirmed)), (np.arange(len(np.log10(xi_nonres_confirmed)))+1.)/np.float(len(np.log10(xi_nonres_confirmed))), drawstyle='steps-post', color='g', ls='--')
#plt.annotate(s='', xy=(P_KS_pos, sum(P_obs_flat <= P_KS_pos)/np.float(len(P_obs_flat))), xytext=(P_KS_pos, sum(P_confirmed <= P_KS_pos)/np.float(len(P_confirmed))), arrowprops=dict(arrowstyle='<->'))
plt.figtext(0.925, 0.35, r'$\mathcal{D} = %s$' % np.round(logxi_KS, 3), ha='right', fontsize=16)
plt.figtext(0.925, 0.3, r'$\mathcal{D} = %s$' % np.round(logxi_res_KS, 3), ha='right', fontsize=16)
plt.figtext(0.925, 0.25, r'$\mathcal{D} = %s$' % np.round(logxi_nonres_KS, 3), ha='right', fontsize=16)
ax.tick_params(axis='both', labelsize=20)
plt.xlim([-0.5, 0.5])
plt.xlabel(r'$\log{\xi}$', fontsize=20)
plt.ylabel('CDF', fontsize=20)
plt.legend(loc='upper left', bbox_to_anchor=(0.01,0.99), ncol=1, frameon=False, fontsize=16) #show the legend
if savefigures == True:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_logxi_CDFs.pdf')
    plt.close()

#To plot the log(xi) CDFs within res:
fig = plt.figure(figsize=(8,4))
plot = GridSpec(1,1,left=0.15,bottom=0.2,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
plt.plot(np.sort(np.log10(xi_res32_obs_flat)), (np.arange(len(np.log10(xi_res32_obs_flat)))+1.)/np.float(len(np.log10(xi_res32_obs_flat))), drawstyle='steps-post', color='r', linewidth=linewidth, label='Near 3:2 MMR')
plt.plot(np.sort(np.log10(xi_res32_confirmed)), (np.arange(len(np.log10(xi_res32_confirmed)))+1.)/np.float(len(np.log10(xi_res32_confirmed))), drawstyle='steps-post', color='r', ls='--')
plt.plot(np.sort(np.log10(xi_res21_obs_flat)), (np.arange(len(np.log10(xi_res21_obs_flat)))+1.)/np.float(len(np.log10(xi_res21_obs_flat))), drawstyle='steps-post', color='b', linewidth=linewidth, label='Near 2:1 MMR')
plt.plot(np.sort(np.log10(xi_res21_confirmed)), (np.arange(len(np.log10(xi_res21_confirmed)))+1.)/np.float(len(np.log10(xi_res21_confirmed))), drawstyle='steps-post', color='b', ls='--')
#plt.annotate(s='', xy=(P_KS_pos, sum(P_obs_flat <= P_KS_pos)/np.float(len(P_obs_flat))), xytext=(P_KS_pos, sum(P_confirmed <= P_KS_pos)/np.float(len(P_confirmed))), arrowprops=dict(arrowstyle='<->'))
plt.figtext(0.925, 0.3, r'$\mathcal{D} = %s$' % np.round(logxi_res32_KS, 3), ha='right', fontsize=16)
plt.figtext(0.925, 0.25, r'$\mathcal{D} = %s$' % np.round(logxi_res21_KS, 3), ha='right', fontsize=16)
ax.tick_params(axis='both', labelsize=20)
plt.xlim([-0.5, 0.5])
plt.xlabel(r'$\log{\xi}$', fontsize=20)
plt.ylabel('CDF', fontsize=20)
plt.legend(loc='upper left', bbox_to_anchor=(0.01,0.99), ncol=1, frameon=False, fontsize=16) #show the legend
if savefigures == True:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_logxi_res_CDFs.pdf')
    plt.close()

plt.show()
plt.close()
#'''





##### To plot the observed multi-systems by period to visualize the systems (similar to Fig 1 in Fabrycky et al. 2014):
#'''
subdirectory = 'Paper_Figures/'

N_multi = sum(Mtot_obs >= 3) #number of simulated multi-systems with 3 or more planets
N_multi_confirmed = sum(M_confirmed >= 3)

i_sorted_P0 = np.argsort(P_obs[Mtot_obs >= 3,0]) #array of indices that would sort the arrays of multi-systems by the innermost period of each system
P_obs_multi = P_obs[Mtot_obs >= 3][i_sorted_P0]
radii_obs_multi = radii_obs[Mtot_obs >= 3][i_sorted_P0]

i_sorted_P0_confirmed = np.argsort(P_sys_confirmed[M_confirmed >= 3,0]) #array of indices that would sort the arrays of multi-systems by the innermost period of each system
P_obs_multi_confirmed = P_sys_confirmed[M_confirmed >= 3][i_sorted_P0_confirmed]
radii_obs_multi_confirmed = radii_sys_confirmed[M_confirmed >= 3][i_sorted_P0_confirmed]

N_sys_per_plot = 170 #number of systems to plot per figure
for i in range(int(np.ceil(float(N_multi)/N_sys_per_plot))):
    fig = plt.figure(figsize=(10,10))
    plot = GridSpec(1,2,left=0.025,bottom=0.1,right=0.975,top=0.95,wspace=0,hspace=0.1)
    
    ax = plt.subplot(plot[0,0])
    plt.title('Kepler multi-planet systems', fontsize=12)
    for j in range(len(P_obs_multi_confirmed[i*N_sys_per_plot:(i+1)*N_sys_per_plot])):
        P_sys = P_obs_multi_confirmed[i*N_sys_per_plot + j]
        radii_sys = radii_obs_multi_confirmed[i*N_sys_per_plot + j]
        P_sys = P_sys[P_sys > 0]
        radii_sys = radii_sys[radii_sys > 0]
        plt.scatter(P_sys, np.ones(len(P_sys))+j, c=np.argsort(radii_sys), s=2.*radii_sys**2.)
        if (j+1)%10 == 0:
            plt.axhline(y=j+1, lw=0.05, color='k')
    plt.gca().set_xscale("log")
    ax.set_xticks([3,10,30,100,300])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_yticks([])
    plt.xlim([2., 500.])
    plt.ylim([0., N_sys_per_plot])
    plt.xlabel(r'$P$ (days)', fontsize=12)

    ax = plt.subplot(plot[0,1])
    plt.title('Simulated multi-planet systems', fontsize=12)
    for j in range(len(P_obs_multi[i*N_sys_per_plot:(i+1)*N_sys_per_plot])):
        P_sys = P_obs_multi[i*N_sys_per_plot + j]
        radii_sys = radii_obs_multi[i*N_sys_per_plot + j]
        P_sys = P_sys[P_sys > 0]
        radii_sys = radii_sys[radii_sys > 0]
        plt.scatter(P_sys, np.ones(len(P_sys))+j, c=np.argsort(radii_sys), s=2.*radii_sys**2.)
        if (j+1)%10 == 0:
            plt.axhline(y=j+1, lw=0.05, color='k')
    plt.gca().set_xscale("log")
    ax.set_xticks([3,10,30,100,300])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_yticks([])
    plt.xlim([2., 500.])
    plt.ylim([0., N_sys_per_plot])
    plt.xlabel(r'$P$ (days)', fontsize=12)

    #if savefigures == True:
        #plt.savefig(savefigures_directory + subdirectory + model_name + '_multis.pdf')
    plt.show()
    plt.close()
#'''







