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

#To define some useful constants and functions:
N_Kep = 150061 #number of Kepler targets satisfying our cuts to give our observed catalog

AU = 1.496*10.**13. #AU in cm
Msun = 1.989*10.**30. #Solar mass in kg
Rsun = 6.957*10.**10. #Solar radius in cm
Mearth = 5.972*10.**24 #Earth mass in kg
Rearth = 6.371*10.**8. #Earth radius in cm

def a_from_P(P, Mstar):
    #This function converts period (days) to semi-major axis (AU) assuming mass of planet m << Mstar (Msun)
    y = (P/365.25)**(2./3.)*(Mstar/1.0)**(1./3.)
    return y

def P_from_a(a, Mstar):
    #This function converts semi-major axis (AU) to period (days) assuming mass of planet m << Mstar (Msun)
    y = 365.25*(a**(3./2.))*(Mstar/1.0)**(-1./2.)
    return y



savefigures = False

loadfiles_directory = 'Clustering_Method_Figures/ExoplanetsSysSim/Non_clustered/Sim6/' #'ExoplanetsSysSim.jl-master/examples/clusters/'; 'Clustering_Method_Figures/ExoplanetsSysSim/Power_law_r1_r2_sigma_r/'
savefigures_directory = 'Clustering_Method_Figures/ExoplanetsSysSim/Non_clustered/'
run_number = ''

model_name = 'ExoplanetsSysSim_Non_clustered_Model' + run_number #'ExoplanetsSysSim_Non_clustered_Model'





##### To load the underlying populations:

#To first read the number of simulated targets and bounds for the periods and radii:
with open(loadfiles_directory + 'periods_all%s.out' % run_number, 'r') as file:
    for line in file:
        if line[:26] == '# num_targets_sim_pass_one':
            N_sim = int(line[28:])
        elif line[:12] == '# min_period':
            P_min = float(line[14:])
        elif line[:12] == '# max_period':
            P_max = float(line[14:])
        elif line[:12] == '# min_radius':
            radii_min = float(line[24:])
        elif line[:12] == '# max_radius':
            radii_max = float(line[24:])

P_per_sys = [] #list to be filled with lists of all periods per system (days)
N_sys_with_planets = 0 #counter for the number of simulated systems with planets
with open(loadfiles_directory + 'periods_all%s.out' % run_number, 'r') as file:
    for line in file:
        if line[0] != '#':
            N_sys_with_planets += 1
            line = line[1:-2].split(', ')
            P_sys = [float(i) for i in line]
            P_per_sys.append(P_sys)
            #print P_sys

e_per_sys = [] #list to be filled with lists of all eccentricities per system
with open(loadfiles_directory + 'eccentricities_all%s.out' % run_number, 'r') as file:
    for line in file:
        if line[0] != '#':
            line = line[1:-2].split(', ')
            e_sys = [float(i) for i in line]
            e_per_sys.append(e_sys)
            #print e_sys

radii_per_sys = [] #list to be filled with lists of all planet radii per system (solar radii)
with open(loadfiles_directory + 'radii_all%s.out' % run_number, 'r') as file:
    for line in file:
        if line[0] != '#':
            line = line[1:-2].split(', ')
            radii_sys = [float(i) for i in line]
            radii_per_sys.append(radii_sys)
            #print radii_sys

mass_per_sys = [] #list to be filled with lists of all planet radii per system (solar masses)
with open(loadfiles_directory + 'masses_all%s.out' % run_number, 'r') as file:
    for line in file:
        if line[0] != '#':
            line = line[1:-2].split(', ')
            mass_sys = [float(i) for i in line]
            mass_per_sys.append(mass_sys)
            #print mass_sys

stellar_mass_all = np.loadtxt(loadfiles_directory + 'stellar_masses_with_planets%s.out' % run_number) #array of stellar masses of all the systems with a planetary system, in solar masses
stellar_radii_all = np.loadtxt(loadfiles_directory + 'stellar_radii_with_planets%s.out' % run_number) #array of stellar radii of all the systems with a planetary system, in solar radii

P_all = [] #list to be zero-padded so each list of periods is sorted and has the same length, and then converted to an array
e_all = [] #list to be zero-padded so each list of eccentricities is sorted (by period) and has the same length, and then converted to an array
radii_all = [] #list to be zero-padded so each list of radii is sorted (by period) and has the same length, and then converted to an array
mass_all = [] #list to be zero-padded so each list of masses is sorted (by period) and has the same length, and then converted to an array

Pmin = 0. #set a minimum period (days), discarding planets less than this period

Mmax = max(len(x) for x in P_per_sys) #maximum planet multiplicity generated by the clustering method
for i in range(len(P_per_sys)):
    i_sorted = np.argsort(P_per_sys[i]) #array of indices which would sort the system by period
    P_sorted = np.array(P_per_sys[i])[i_sorted]
    P_sorted_cut = P_sorted[P_sorted > Pmin]
    e_sorted_cut = np.array(e_per_sys[i])[i_sorted][P_sorted > Pmin]
    radii_sorted_cut = np.array(radii_per_sys[i])[i_sorted][P_sorted > Pmin]
    mass_sorted_cut = np.array(mass_per_sys[i])[i_sorted][P_sorted > Pmin]
    
    P_sys = list(P_sorted_cut) + [0]*(Mmax - len(P_sorted_cut)) #zero-pad the list up to Mmax elements
    e_sys = list(e_sorted_cut) + [0]*(Mmax - len(e_sorted_cut)) #zero-pad the list up to Mmax elements
    radii_sys = list(radii_sorted_cut) + [0]*(Mmax - len(radii_sorted_cut)) #zero-pad the list up to Mmax elements
    mass_sys = list(mass_sorted_cut) + [0]*(Mmax - len(mass_sorted_cut)) #zero-pad the list up to Mmax elements

    P_all.append(P_sys)
    e_all.append(e_sys)
    radii_all.append(radii_sys)
    mass_all.append(mass_sys)
P_all = np.array(P_all)
e_all = np.array(e_all)
radii_all = np.array(radii_all)
mass_all = np.array(mass_all)

#To convert the radii and masses to Earth units:
radii_all = radii_all*(Rsun/Rearth) #radii in Earth radii
mass_all = mass_all*(Msun/Mearth) #masses in Earth masses

Mtot_all = np.sum(P_all > 0, axis=1) #array of all planet multiplicites





#To calculate the underlying period ratios, radii ratios, and separations in mutual Hill radii:
Rm_all = [] #list to be filled with all the period ratios
radii_ratio_all = [] #list to be filled with all the radii ratios
N_mH_all = [] #list to be filled with all the separations between adjacent planet pairs in units of mutual Hill radii
for i in range(len(P_all)):
    Mstar_system = stellar_mass_all[i] #mass of the star for this system, in solar masses
    P_all_system = P_all[i][P_all[i] > 0]
    e_all_system = e_all[i][P_all[i] > 0]
    radii_all_system = radii_all[i][P_all[i] > 0]
    mass_all_system = mass_all[i][P_all[i] > 0]
    
    #To calculate all the period ratios:
    Rm_all_system = list(P_all_system[1:]/P_all_system[0:-1]) #list of period ratios in this system
    Rm_all_system = np.array(Rm_all_system + [0]*(Mmax - 1 - len(Rm_all_system))) #to add filler 0's to Rm_all_system to pad it to Mmax - 1 elements
    Rm_all.append(Rm_all_system)

    #To calculate all the radii ratios:
    radii_ratio_all_system = list(radii_all_system[1:]/radii_all_system[0:-1]) #list of radii ratios in this system
    radii_ratio_all_system = np.array(radii_ratio_all_system + [0]*(Mmax - 1 - len(radii_ratio_all_system))) #to add filler 0's to radii_ratio_all_system to pad it to Mmax - 1 elements
    radii_ratio_all.append(radii_ratio_all_system)

    #To calculate all the separations in mutual Hill radii between adjacent planet pairs:
    a_all_system = a_from_P(P_all_system, Mstar_system)
    R_mH_all_system = ((a_all_system[0:-1] + a_all_system[1:])/2.)*(Mearth*(mass_all_system[0:-1] + mass_all_system[1:])/(3.*Mstar_system*Msun))**(1./3.) #mutual Hill radii between adjacent planet pairs in this system, in AU
    R_sep_all_system = a_all_system[1:] - a_all_system[0:-1] #separations between adjacent planet pairs in this system, in AU, ignoring eccentricities
    N_mH_all_system = list(R_sep_all_system/R_mH_all_system) #separations between adjacent planet pairs in this system, in mutual Hill radii
    N_mH_all_system = np.array(N_mH_all_system + [0]*(Mmax - 1 - len(N_mH_all_system))) #to add filler 0's to N_mH_all_system to pad it to Mmax - 1 elements
    N_mH_all.append(N_mH_all_system)


Rm_all = np.array(Rm_all)
radii_ratio_all = np.array(radii_ratio_all)
N_mH_all = np.array(N_mH_all)

P_all_flat = P_all.flatten() #all the periods of all the planets
P_all_flat = P_all_flat[P_all_flat > 0]

Rm_all_flat = Rm_all.flatten() #all the period ratios of all the adjacent planets
Rm_all_flat = Rm_all_flat[Rm_all_flat > 0]

N_mH_all_flat = N_mH_all.flatten() #all the separations of all the adjacent planets in units of their mutual Hill radii
N_mH_all_flat = N_mH_all_flat[N_mH_all_flat > 0]

e_all_flat = e_all.flatten() #all the eccentricities of all the planets
e_all_flat = e_all_flat[e_all_flat > 0]

radii_all_flat = radii_all.flatten() #all the planet radii, in Earth radii
radii_all_flat = radii_all_flat[radii_all_flat > 0]

radii_ratio_all_flat = radii_ratio_all.flatten() #all the radii ratios of all the adjacent planets
radii_ratio_all_flat = radii_ratio_all_flat[radii_ratio_all_flat > 0]

mass_all_flat = mass_all.flatten() #all the planet masses, in Earth masses
mass_all_flat = mass_all_flat[mass_all_flat > 0]





#To plot the underlying distributions (total planet multiplicity, periods, period ratios):
'''
fig = plt.figure(figsize=(16,8))
plot = GridSpec(2,2,left=0.075,bottom=0.115,right=0.95,top=0.925,wspace=0.15,hspace=0.4)

ax = plt.subplot(plot[0,0])
plt.title('Underlying planet multiplicity', fontsize=24)
x = Mtot_all[Mtot_all > 0]
counts, bins = np.histogram(x, bins=np.max(x)+1, range=(-0.5, np.max(x)+0.5))
bins_mid = (bins[:-1] + bins[1:])/2.
plt.plot(bins_mid, counts, 'o-', color='k', label='%s systems with planets' % len(x))
ax.tick_params(axis='both', labelsize=20)
plt.xlim([np.min(x), np.max(x)])
plt.xlabel(r'$M_{\rm tot}$', fontsize=20)
plt.ylabel('Number', fontsize=20)
plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[0,1])
plt.title('All periods', fontsize=24)
hist = plt.hist(P_all_flat, bins=np.logspace(np.log10(np.min(P_all_flat)), np.log10(np.max(P_all_flat)), 101), histtype='step', log=True, color='k', label=r'All planets')
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([np.min(P_all_flat), 1.1*np.max(P_all_flat)])
plt.xlabel(r'$P$ (days)', fontsize=20)
#plt.ylabel('Number', fontsize=20)
#plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[1,0])
plt.title('All period ratios', fontsize=24)
x = Rm_all_flat
plt.hist(x, bins=np.logspace(np.log10(1.), np.log10(np.max(x)), 101), histtype='step', color='k', label='All adjacent pairs')
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
#plt.xlim([1,5])
plt.xlabel(r'$P_{i+1}/P_i$', fontsize=20)
plt.ylabel('Number', fontsize=20)
#plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[1,1])
plt.title('All period ratios ($P_{i+1}/P_i < 5$)', fontsize=24)
x = Rm_all_flat[Rm_all_flat < 5]
plt.hist(x, bins=np.logspace(np.log10(1.), np.log10(5.), 101), histtype='step', color='k', label='All adjacent pairs')
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([1,5])
plt.xlabel(r'$P_{i+1}/P_i$', fontsize=20)
#plt.ylabel('Number', fontsize=20)
#plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

if savefigures == True:
    plt.savefig(savefigures_directory + model_name + '_underlying_periods.pdf')
else:
    plt.show()
plt.close()
'''

#To plot the underlying distributions (eccentricities, masses, radii, radii ratios):
'''
fig = plt.figure(figsize=(16,8))
plot = GridSpec(2,2,left=0.075,bottom=0.115,right=0.95,top=0.925,wspace=0.15,hspace=0.4)

ax = plt.subplot(plot[0,0])
plt.title('Orbital eccentricities', fontsize=24)
x = e_all_flat
plt.hist(x, bins=np.linspace(0,0.2,101), histtype='step', color='k', label=r'All planets')
ax.tick_params(axis='both', labelsize=20)
plt.xlim([np.min(x), np.max(x)])
plt.xlabel(r'$e$', fontsize=20)
plt.ylabel('Number', fontsize=20)
plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[0,1])
plt.title('Planet masses', fontsize=24)
x = mass_all_flat
plt.hist(x, bins=np.logspace(np.log10(np.min(x)),np.log10(np.max(x)),101), histtype='step', color='k', label=r'All planets')
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([np.min(mass_all_flat), 1.1*np.max(mass_all_flat)])
plt.xlabel(r'$M_p$ ($M_\oplus$)', fontsize=20)
#plt.ylabel('Number', fontsize=20)
#plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[1,0])
plt.title('Planet radii', fontsize=24)
x = radii_all_flat
plt.hist(x, bins=np.logspace(np.log10(np.min(x)),np.log10(np.max(x)),101), histtype='step', color='k', label='All planets')
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([0.5,10])
plt.xlabel(r'$R_p$ ($R_\oplus$)', fontsize=20)
plt.ylabel('Number', fontsize=20)
#plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

ax = plt.subplot(plot[1,1])
plt.title('Planet radii ratios', fontsize=24)
x = radii_ratio_all_flat
plt.hist(x, bins=np.logspace(np.log10(np.min(x)),np.log10(np.max(x)),101), histtype='step', color='k', label='All adjacent pairs')
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
#plt.xlim([np.min(x),np.max(x)])
plt.xlabel(r'$R_{p,i+1}/R_{p,i}$', fontsize=20)
#plt.ylabel('Number', fontsize=20)
#plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

if savefigures == True:
    plt.savefig(savefigures_directory + model_name + '_underlying_properties.pdf')
else:
    plt.show()
plt.close()
'''

#To plot the underlying distributions (separations in mutual Hill radii):
'''
fig = plt.figure(figsize=(16,8))
plot = GridSpec(2,2,left=0.075,bottom=0.115,right=0.95,top=0.925,wspace=0.15,hspace=0.4)

ax = plt.subplot(plot[0,0])
plt.title('Separation in mutual Hill radii', fontsize=24)
x = N_mH_all_flat
plt.hist(x, bins=np.logspace(np.log10(np.min(x)),np.log10(np.max(x)),101), histtype='step', color='k', label=r'All planets')
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
#plt.xlim([np.min(x), np.max(x)])
plt.xlabel(r'$\Delta$', fontsize=20)
plt.ylabel('Number', fontsize=20)
plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=16) #show the legend

if savefigures == True:
    plt.savefig(savefigures_directory + model_name + '_underlying_stability.pdf')
else:
    plt.show()
plt.close()
'''





##### To plot the underlying multi-systems by period to visualize the systems (similar to Fig 1 in Fabrycky et al. 2014, but for ALL the planets):
##### Note: since there are way too many simulated systems to plot them all, we will randomly sample a number of systems to plot
'''
N_multi = sum(Mtot_all >= 3) #number of simulated multi-systems with 3 or more planets
i_multi = np.arange(len(Mtot_all))[Mtot_all >= 3] #array of indices of all multi-systems with 3 or more planets

N_sys_per_plot = 100 #number of systems to sample and plot per figure
i_multi_sample = np.random.choice(i_multi, N_sys_per_plot, replace=False) #array of indices of a sample of multi-systems with 3 or more planets

i_sorted_P0 = np.argsort(P_all[i_multi_sample,0]) #array of indices that would sort the arrays of the sample of multi-systems by the innermost period of each system
P_sample_multi = P_all[i_multi_sample][i_sorted_P0]
radii_sample_multi = radii_all[i_multi_sample][i_sorted_P0]

fig = plt.figure(figsize=(10,10))
plot = GridSpec(1,1,left=0.05,bottom=0.2,right=0.95,top=0.95,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
for j in range(len(P_sample_multi)):
    P_sys = P_sample_multi[j]
    radii_sys = radii_sample_multi[j]
    P_sys = P_sys[P_sys > 0]
    radii_sys = radii_sys[radii_sys > 0]
    plt.scatter(P_sys, np.ones(len(P_sys))+j, c=np.argsort(radii_sys), s=10.*radii_sys**2.)
    plt.axhline(y=j+1, lw=0.1, color='k')
plt.gca().set_xscale("log")
ax.set_yticks([])
plt.xlim([2., 500.])
plt.ylim([0., N_sys_per_plot])
plt.xlabel(r'$P$ (days)', fontsize=20)
plt.show()
plt.close()
'''










##### If we want to load and simultaneously plot another model, for model comparison:

loadfiles_directory2 = 'Clustering_Method_Figures/ExoplanetsSysSim/Power_law_r1_r2_sigma_r/Sim2/' #'ExoplanetsSysSim.jl-master/examples/clusters/'; 'Clustering_Method_Figures/ExoplanetsSysSim/Power_law_r1_r2_sigma_r/'
savefigures_directory2 = 'Clustering_Method_Figures/ExoplanetsSysSim/Power_law_r1_r2_sigma_r/'
run_number2 = ''

#THIS REWRITES THE MODEL NAME:
model_name = 'ExoplanetsSysSim_Models_compare' #'ExoplanetsSysSim_Non_clustered_Model'





##### To load the underlying populations:

#To first read the number of simulated targets and bounds for the periods and radii:
with open(loadfiles_directory2 + 'periods_all%s.out' % run_number2, 'r') as file:
    for line in file:
        if line[:26] == '# num_targets_sim_pass_one':
            N_sim2 = int(line[28:])
        ##### Assuming the bounds on the periods and radii are the same as the previously loaded model!

P_per_sys2 = [] #list to be filled with lists of all periods per system (days)
N_sys_with_planets2 = 0 #counter for the number of simulated systems with planets
with open(loadfiles_directory2 + 'periods_all%s.out' % run_number2, 'r') as file:
    for line in file:
        if line[0] != '#':
            N_sys_with_planets2 += 1
            line = line[1:-2].split(', ')
            P_sys = [float(i) for i in line]
            P_per_sys2.append(P_sys)

e_per_sys2 = [] #list to be filled with lists of all eccentricities per system
with open(loadfiles_directory2 + 'eccentricities_all%s.out' % run_number2, 'r') as file:
    for line in file:
        if line[0] != '#':
            line = line[1:-2].split(', ')
            e_sys = [float(i) for i in line]
            e_per_sys2.append(e_sys)

radii_per_sys2 = [] #list to be filled with lists of all planet radii per system (solar radii)
with open(loadfiles_directory2 + 'radii_all%s.out' % run_number2, 'r') as file:
    for line in file:
        if line[0] != '#':
            line = line[1:-2].split(', ')
            radii_sys = [float(i) for i in line]
            radii_per_sys2.append(radii_sys)

mass_per_sys2 = [] #list to be filled with lists of all planet radii per system (solar masses)
with open(loadfiles_directory2 + 'masses_all%s.out' % run_number2, 'r') as file:
    for line in file:
        if line[0] != '#':
            line = line[1:-2].split(', ')
            mass_sys = [float(i) for i in line]
            mass_per_sys2.append(mass_sys)

stellar_mass_all2 = np.loadtxt(loadfiles_directory2 + 'stellar_masses_with_planets%s.out' % run_number2) #array of stellar masses of all the systems with a planetary system, in solar masses
stellar_radii_all2 = np.loadtxt(loadfiles_directory2 + 'stellar_radii_with_planets%s.out' % run_number2) #array of stellar radii of all the systems with a planetary system, in solar radii

P_all2 = [] #list to be zero-padded so each list of periods is sorted and has the same length, and then converted to an array
e_all2 = [] #list to be zero-padded so each list of eccentricities is sorted (by period) and has the same length, and then converted to an array
radii_all2 = [] #list to be zero-padded so each list of radii is sorted (by period) and has the same length, and then converted to an array
mass_all2 = [] #list to be zero-padded so each list of masses is sorted (by period) and has the same length, and then converted to an array

Pmin2 = 0. #set a minimum period (days), discarding planets less than this period

Mmax2 = max(len(x) for x in P_per_sys2) #maximum planet multiplicity generated by the clustering method
for i in range(len(P_per_sys2)):
    i_sorted = np.argsort(P_per_sys2[i]) #array of indices which would sort the system by period
    P_sorted = np.array(P_per_sys2[i])[i_sorted]
    P_sorted_cut = P_sorted[P_sorted > Pmin]
    e_sorted_cut = np.array(e_per_sys2[i])[i_sorted][P_sorted > Pmin]
    radii_sorted_cut = np.array(radii_per_sys2[i])[i_sorted][P_sorted > Pmin]
    mass_sorted_cut = np.array(mass_per_sys2[i])[i_sorted][P_sorted > Pmin]
    
    P_sys = list(P_sorted_cut) + [0]*(Mmax2 - len(P_sorted_cut)) #zero-pad the list up to Mmax elements
    e_sys = list(e_sorted_cut) + [0]*(Mmax2 - len(e_sorted_cut)) #zero-pad the list up to Mmax elements
    radii_sys = list(radii_sorted_cut) + [0]*(Mmax2 - len(radii_sorted_cut)) #zero-pad the list up to Mmax elements
    mass_sys = list(mass_sorted_cut) + [0]*(Mmax2 - len(mass_sorted_cut)) #zero-pad the list up to Mmax elements
    
    P_all2.append(P_sys)
    e_all2.append(e_sys)
    radii_all2.append(radii_sys)
    mass_all2.append(mass_sys)
P_all2 = np.array(P_all2)
e_all2 = np.array(e_all2)
radii_all2 = np.array(radii_all2)
mass_all2 = np.array(mass_all2)

#To convert the radii and masses to Earth units:
radii_all2 = radii_all2*(Rsun/Rearth) #radii in Earth radii
mass_all2 = mass_all2*(Msun/Mearth) #masses in Earth masses

Mtot_all2 = np.sum(P_all2 > 0, axis=1) #array of all planet multiplicites





#To calculate the underlying period ratios, radii ratios, and separations in mutual Hill radii:
Rm_all2 = [] #list to be filled with all the period ratios
radii_ratio_all2 = [] #list to be filled with all the radii ratios
N_mH_all2 = [] #list to be filled with all the separations between adjacent planet pairs in units of mutual Hill radii
for i in range(len(P_all2)):
    Mstar_system = stellar_mass_all2[i] #mass of the star for this system, in solar masses
    P_all_system = P_all2[i][P_all2[i] > 0]
    e_all_system = e_all2[i][P_all2[i] > 0]
    radii_all_system = radii_all2[i][P_all2[i] > 0]
    mass_all_system = mass_all2[i][P_all2[i] > 0]
    
    #To calculate all the period ratios:
    Rm_all_system = list(P_all_system[1:]/P_all_system[0:-1]) #list of period ratios in this system
    Rm_all_system = np.array(Rm_all_system + [0]*(Mmax2 - 1 - len(Rm_all_system))) #to add filler 0's to Rm_all_system to pad it to Mmax - 1 elements
    Rm_all2.append(Rm_all_system)
    
    #To calculate all the radii ratios:
    radii_ratio_all_system = list(radii_all_system[1:]/radii_all_system[0:-1]) #list of radii ratios in this system
    radii_ratio_all_system = np.array(radii_ratio_all_system + [0]*(Mmax2 - 1 - len(radii_ratio_all_system))) #to add filler 0's to radii_ratio_all_system to pad it to Mmax - 1 elements
    radii_ratio_all2.append(radii_ratio_all_system)
    
    #To calculate all the separations in mutual Hill radii between adjacent planet pairs:
    a_all_system = a_from_P(P_all_system, Mstar_system)
    R_mH_all_system = ((a_all_system[0:-1] + a_all_system[1:])/2.)*(Mearth*(mass_all_system[0:-1] + mass_all_system[1:])/(3.*Mstar_system*Msun))**(1./3.) #mutual Hill radii between adjacent planet pairs in this system, in AU
    R_sep_all_system = a_all_system[1:] - a_all_system[0:-1] #separations between adjacent planet pairs in this system, in AU, ignoring eccentricities
    N_mH_all_system = list(R_sep_all_system/R_mH_all_system) #separations between adjacent planet pairs in this system, in mutual Hill radii
    N_mH_all_system = np.array(N_mH_all_system + [0]*(Mmax2 - 1 - len(N_mH_all_system))) #to add filler 0's to N_mH_all_system to pad it to Mmax - 1 elements
    N_mH_all2.append(N_mH_all_system)


Rm_all2 = np.array(Rm_all2)
radii_ratio_all2 = np.array(radii_ratio_all2)
N_mH_all2 = np.array(N_mH_all2)

P_all_flat2 = P_all2.flatten() #all the periods of all the planets
P_all_flat2 = P_all_flat2[P_all_flat2 > 0]

Rm_all_flat2 = Rm_all2.flatten() #all the period ratios of all the adjacent planets
Rm_all_flat2 = Rm_all_flat2[Rm_all_flat2 > 0]

N_mH_all_flat2 = N_mH_all2.flatten() #all the separations of all the adjacent planets in units of their mutual Hill radii
N_mH_all_flat2 = N_mH_all_flat2[N_mH_all_flat2 > 0]

e_all_flat2 = e_all2.flatten() #all the eccentricities of all the planets
e_all_flat2 = e_all_flat2[e_all_flat2 > 0]

radii_all_flat2 = radii_all2.flatten() #all the planet radii, in Earth radii
radii_all_flat2 = radii_all_flat2[radii_all_flat2 > 0]

radii_ratio_all_flat2 = radii_ratio_all2.flatten() #all the radii ratios of all the adjacent planets
radii_ratio_all_flat2 = radii_ratio_all_flat2[radii_ratio_all_flat2 > 0]

mass_all_flat2 = mass_all2.flatten() #all the planet masses, in Earth masses
mass_all_flat2 = mass_all_flat2[mass_all_flat2 > 0]




#'''
#To plot the underlying distributions of the two models as individual panels:

savefigures_directory_compare = 'Clustering_Method_Figures/ExoplanetsSysSim/Underlying_Compare/'
subdirectory = 'Talk_Figures/'

#Make sure these labels match the models being loaded!
model1_label = 'Non clustered'
model2_label = 'Clustered'

#To plot the multiplicities:
fig = plt.figure(figsize=(8,4))
plot = GridSpec(1,1,left=0.15,bottom=0.2,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
x1 = Mtot_all
x2 = Mtot_all2
max_M = np.max((np.max(x1), np.max(x2)))
counts, bins = np.histogram(x1, bins=max_M+1, range=(-0.5, max_M+0.5))
counts[0] = N_sim - len(Mtot_all) #to compute the number of systems with no planets
bins_mid = (bins[:-1] + bins[1:])/2.
plt.plot(bins_mid, counts/float(N_sim), 'o-', color='r', label= model1_label)
counts, bins = np.histogram(x2, bins=max_M+1, range=(-0.5, max_M+0.5))
counts[0] = N_sim2 - len(Mtot_all2) #to compute the number of systems with no planets
plt.plot(bins_mid, counts/float(N_sim2), 'o-', color='b', label= model2_label)
ax.tick_params(axis='both', labelsize=20)
plt.xlim([0, max_M])
plt.xlabel(r'Number of planets', fontsize=20)
plt.ylabel('Fraction', fontsize=20)
plt.legend(loc='upper right', bbox_to_anchor=(0.99,0.99), ncol=1, frameon=False, fontsize=16) #show the legend
if savefigures == True:
    plt.savefig(savefigures_directory_compare + subdirectory + model_name + '_underlying_multiplicities_compare.pdf')
    plt.close()

#To plot the periods:
fig = plt.figure(figsize=(8,4))
plot = GridSpec(1,1,left=0.15,bottom=0.2,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
x1 = P_all_flat
x2 = P_all_flat2
xmin, xmax = np.min((np.min(x1), np.min(x2))), np.max((np.max(x1), np.max(x2)))
hist1 = plt.hist(x1, bins=np.logspace(np.log10(xmin), np.log10(xmax), 101), histtype='step', weights=np.ones(len(x1))/len(x1), log=True, color='r', label=model1_label)
hist2 = plt.hist(x2, bins=np.logspace(np.log10(xmin), np.log10(xmax), 101), histtype='step', weights=np.ones(len(x2))/len(x2), log=True, color='b', label=model2_label)
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([np.min(P_all_flat), 1.1*np.max(P_all_flat)])
plt.ylim([np.min((np.min(hist1[0][hist1[0] > 0]), np.min(hist2[0][hist2[0] > 0]))), np.max((np.max(hist1[0][hist1[0] > 0]), np.max(hist2[0][hist2[0] > 0])))])
plt.xlabel(r'$P$ (days)', fontsize=20)
plt.ylabel('Fraction', fontsize=20)
plt.legend(loc='upper right', bbox_to_anchor=(0.99,0.99), ncol=1, frameon=False, fontsize=16) #show the legend
if savefigures == True:
    plt.savefig(savefigures_directory_compare + subdirectory + model_name + '_underlying_periods_compare.pdf')
    plt.close()

#To plot the period ratios:
fig = plt.figure(figsize=(8,4))
plot = GridSpec(1,1,left=0.15,bottom=0.2,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
x1 = Rm_all_flat
x2 = Rm_all_flat2
xmin, xmax = np.min((np.min(x1), np.min(x2))), np.max((np.max(x1), np.max(x2)))
hist1 = plt.hist(x1, bins=np.logspace(np.log10(xmin), np.log10(xmax), 101), histtype='step', weights=np.ones(len(x1))/len(x1), log=True, color='r', label=model1_label)
hist2 = plt.hist(x2, bins=np.logspace(np.log10(xmin), np.log10(xmax), 101), histtype='step', weights=np.ones(len(x2))/len(x2), log=True, color='b', label=model2_label)
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
#plt.xlim([1,5])
plt.ylim([np.min((np.min(hist1[0][hist1[0] > 0]), np.min(hist2[0][hist2[0] > 0]))), np.max((np.max(hist1[0][hist1[0] > 0]), np.max(hist2[0][hist2[0] > 0])))])
plt.xlabel(r'$P_{i+1}/P_i$', fontsize=20)
plt.ylabel('Fraction', fontsize=20)
plt.legend(loc='upper right', bbox_to_anchor=(0.99,0.99), ncol=1, frameon=False, fontsize=16) #show the legend
if savefigures == True:
    plt.savefig(savefigures_directory_compare + subdirectory + model_name + '_underlying_periodratios_compare.pdf')
    plt.close()

#To plot the eccentricities:
fig = plt.figure(figsize=(8,4))
plot = GridSpec(1,1,left=0.15,bottom=0.2,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
x1 = e_all_flat
x2 = e_all_flat2
xmin, xmax = np.min((np.min(x1), np.min(x2))), np.max((np.max(x1), np.max(x2)))
hist1 = plt.hist(x1, bins=np.linspace(xmin,xmax,101), histtype='step', weights=np.ones(len(x1))/len(x1), color='r', label=model1_label)
hist2 = plt.hist(x2, bins=np.linspace(xmin,xmax,101), histtype='step', weights=np.ones(len(x2))/len(x2), color='b', label=model2_label)
ax.tick_params(axis='both', labelsize=20)
plt.ylim([np.min((np.min(hist1[0][hist1[0] > 0]), np.min(hist2[0][hist2[0] > 0]))), np.max((np.max(hist1[0][hist1[0] > 0]), np.max(hist2[0][hist2[0] > 0])))])
plt.xlabel(r'$e$', fontsize=20)
plt.ylabel('Fraction', fontsize=20)
plt.legend(loc='upper right', bbox_to_anchor=(0.99,0.99), ncol=1, frameon=False, fontsize=16) #show the legend
if savefigures == True:
    plt.savefig(savefigures_directory_compare + subdirectory + model_name + '_underlying_eccentricities_compare.pdf')
    plt.close()

#To plot the masses:
fig = plt.figure(figsize=(8,4))
plot = GridSpec(1,1,left=0.15,bottom=0.2,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
x1 = mass_all_flat
x2 = mass_all_flat2
xmin, xmax = np.min((np.min(x1), np.min(x2))), np.max((np.max(x1), np.max(x2)))
hist1 = plt.hist(x1, bins=np.logspace(np.log10(xmin), np.log10(xmax), 101), histtype='step', weights=np.ones(len(x1))/len(x1), log=True, color='r', label=model1_label)
hist2 = plt.hist(x2, bins=np.logspace(np.log10(xmin), np.log10(xmax), 101), histtype='step', weights=np.ones(len(x2))/len(x2), log=True, color='b', label=model2_label)
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
#plt.xlim([np.min(mass_all_flat), 1.1*np.max(mass_all_flat)])
plt.ylim([np.min((np.min(hist1[0][hist1[0] > 0]), np.min(hist2[0][hist2[0] > 0]))), np.max((np.max(hist1[0][hist1[0] > 0]), np.max(hist2[0][hist2[0] > 0])))])
plt.xlabel(r'$M_p$ ($M_\oplus$)', fontsize=20)
plt.ylabel('Fraction', fontsize=20)
plt.legend(loc='upper right', bbox_to_anchor=(0.99,0.99), ncol=1, frameon=False, fontsize=16) #show the legend
if savefigures == True:
    plt.savefig(savefigures_directory_compare + subdirectory + model_name + '_underlying_masses_compare.pdf')
    plt.close()

#To plot the radii:
fig = plt.figure(figsize=(8,4))
plot = GridSpec(1,1,left=0.15,bottom=0.2,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
x1 = radii_all_flat
x2 = radii_all_flat2
xmin, xmax = np.min((np.min(x1), np.min(x2))), np.max((np.max(x1), np.max(x2)))
hist1 = plt.hist(x1, bins=np.logspace(np.log10(xmin), np.log10(xmax), 101), histtype='step', weights=np.ones(len(x1))/len(x1), log=True, color='r', label=model1_label)
hist2 = plt.hist(x2, bins=np.logspace(np.log10(xmin), np.log10(xmax), 101), histtype='step', weights=np.ones(len(x2))/len(x2), log=True, color='b', label=model2_label)
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([0.5,10])
plt.ylim([np.min((np.min(hist1[0][hist1[0] > 0]), np.min(hist2[0][hist2[0] > 0]))), np.max((np.max(hist1[0][hist1[0] > 0]), np.max(hist2[0][hist2[0] > 0])))])
plt.xlabel(r'$R_p$ ($R_\oplus$)', fontsize=20)
plt.ylabel('Fraction', fontsize=20)
plt.legend(loc='upper right', bbox_to_anchor=(0.99,0.99), ncol=1, frameon=False, fontsize=16) #show the legend
if savefigures == True:
    plt.savefig(savefigures_directory_compare + subdirectory + model_name + '_underlying_radii_compare.pdf')
    plt.close()

#To plot the radii ratios:
fig = plt.figure(figsize=(8,4))
plot = GridSpec(1,1,left=0.15,bottom=0.2,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
x1 = radii_ratio_all_flat
x2 = radii_ratio_all_flat2
xmin, xmax = np.min((np.min(x1), np.min(x2))), np.max((np.max(x1), np.max(x2)))
hist1 = plt.hist(x1, bins=np.logspace(np.log10(xmin), np.log10(xmax), 101), histtype='step', weights=np.ones(len(x1))/len(x1), log=False, color='r', label=model1_label)
hist2 = plt.hist(x2, bins=np.logspace(np.log10(xmin), np.log10(xmax), 101), histtype='step', weights=np.ones(len(x2))/len(x2), log=False, color='b', label=model2_label)
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.xlim([np.min((np.min(x1), np.min(x2))), np.max((np.max(x1), np.max(x2)))])
plt.ylim([np.min((np.min(hist1[0][hist1[0] > 0]), np.min(hist2[0][hist2[0] > 0]))), np.max((np.max(hist1[0][hist1[0] > 0]), np.max(hist2[0][hist2[0] > 0])))])
plt.xlabel(r'$R_{p,i+1}/R_{p,i}$', fontsize=20)
plt.ylabel('Fraction', fontsize=20)
plt.legend(loc='upper right', bbox_to_anchor=(0.99,0.99), ncol=1, frameon=False, fontsize=16) #show the legend
if savefigures == True:
    plt.savefig(savefigures_directory_compare + subdirectory + model_name + '_underlying_radiiratios_compare.pdf')
    plt.close()

#To plot the separations in mutual Hill radii:
fig = plt.figure(figsize=(8,4))
plot = GridSpec(1,1,left=0.15,bottom=0.2,right=0.95,top=0.925,wspace=0.1,hspace=0.1)
ax = plt.subplot(plot[0,0])
x1 = N_mH_all_flat
x2 = N_mH_all_flat2
xmin, xmax = np.min((np.min(x1), np.min(x2))), np.max((np.max(x1), np.max(x2)))
hist1 = plt.hist(x1, bins=np.logspace(np.log10(xmin), np.log10(xmax), 101), histtype='step', weights=np.ones(len(x1))/len(x1), log=True, color='r', label=model1_label)
hist2 = plt.hist(x2, bins=np.logspace(np.log10(xmin), np.log10(xmax), 101), histtype='step', weights=np.ones(len(x2))/len(x2), log=True, color='b', label=model2_label)
plt.gca().set_xscale("log")
ax.tick_params(axis='both', labelsize=20)
plt.ylim([np.min((np.min(hist1[0][hist1[0] > 0]), np.min(hist2[0][hist2[0] > 0]))), np.max((np.max(hist1[0][hist1[0] > 0]), np.max(hist2[0][hist2[0] > 0])))])
plt.xlabel(r'$\Delta$', fontsize=20)
plt.ylabel('Fraction', fontsize=20)
plt.legend(loc='upper right', bbox_to_anchor=(0.99,0.99), ncol=1, frameon=False, fontsize=16) #show the legend
if savefigures == True:
    plt.savefig(savefigures_directory_compare + subdirectory + model_name + '_underlying_stability_compare.pdf')
    plt.close()

plt.show()
plt.close()
#'''





##### To plot the underlying multi-systems by period to visualize the systems (similar to Fig 1 in Fabrycky et al. 2014, but for ALL the planets):
##### Note: since there are way too many simulated systems to plot them all, we will randomly sample a number of systems to plot
'''
N_multi = sum(Mtot_all >= 3) #number of simulated multi-systems with 3 or more planets
N_multi2 = sum(Mtot_all2 >= 3) #number of simulated multi-systems with 3 or more planets
i_multi = np.arange(len(Mtot_all))[Mtot_all >= 3] #array of indices of all multi-systems with 3 or more planets
i_multi2 = np.arange(len(Mtot_all2))[Mtot_all2 >= 3] #array of indices of all multi-systems with 3 or more planets

N_sys_per_plot = 100 #number of systems to sample and plot per figure
i_multi_sample = np.random.choice(i_multi, N_sys_per_plot, replace=False) #array of indices of a sample of multi-systems with 3 or more planets
i_multi_sample2 = np.random.choice(i_multi2, N_sys_per_plot, replace=False) #array of indices of a sample of multi-systems with 3 or more planets

i_sorted_P0 = np.argsort(P_all[i_multi_sample,0]) #array of indices that would sort the arrays of the sample of multi-systems by the innermost period of each system
i_sorted_P0_2 = np.argsort(P_all2[i_multi_sample2,0]) #array of indices that would sort the arrays of the sample of multi-systems by the innermost period of each system

P_sample_multi = P_all[i_multi_sample][i_sorted_P0]
radii_sample_multi = radii_all[i_multi_sample][i_sorted_P0]
P_sample_multi2 = P_all2[i_multi_sample2][i_sorted_P0_2]
radii_sample_multi2 = radii_all2[i_multi_sample2][i_sorted_P0_2]

fig = plt.figure(figsize=(10,10))
plot = GridSpec(1,2,left=0.05,bottom=0.2,right=0.95,top=0.95,wspace=0,hspace=0)

ax = plt.subplot(plot[0,0])
plt.title(model1_label, fontsize=20)
for j in range(len(P_sample_multi)):
    P_sys = P_sample_multi[j]
    radii_sys = radii_sample_multi[j]
    P_sys = P_sys[P_sys > 0]
    radii_sys = radii_sys[radii_sys > 0]
    plt.scatter(P_sys, np.ones(len(P_sys))+j, c=np.argsort(radii_sys), s=10.*radii_sys**2.)
    plt.axhline(y=j+1, lw=0.1, color='k')
plt.gca().set_xscale("log")
ax.set_yticks([])
plt.xlim([2., 500.])
plt.ylim([0., N_sys_per_plot])
plt.xlabel(r'$P$ (days)', fontsize=20)

ax = plt.subplot(plot[0,1])
plt.title(model2_label, fontsize=20)
for j in range(len(P_sample_multi2)):
    P_sys = P_sample_multi2[j]
    radii_sys = radii_sample_multi2[j]
    P_sys = P_sys[P_sys > 0]
    radii_sys = radii_sys[radii_sys > 0]
    plt.scatter(P_sys, np.ones(len(P_sys))+j, c=np.argsort(radii_sys), s=10.*radii_sys**2.)
    plt.axhline(y=j+1, lw=0.1, color='k')
plt.gca().set_xscale("log")
ax.set_yticks([])
plt.xlim([2., 500.])
plt.ylim([0., N_sys_per_plot])
plt.xlabel(r'$P$ (days)', fontsize=20)

plt.show()
plt.close()
'''

