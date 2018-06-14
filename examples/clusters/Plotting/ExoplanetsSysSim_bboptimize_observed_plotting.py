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





##### This module will be used to plot simulated catalogs generated from ExoplanetsSysSim using the best active parameter values from optimization runs:

#To define some useful constants:
N_Kep = 150061 #number of Kepler targets satisfying our cuts to give our observed catalog

AU = 1.496*10.**13. #AU in cm
Msun = 1.989*10.**30. #Solar mass in kg
Rsun = 6.957*10.**10. #Solar radius in cm
Rearth = 6.371*10.**8. #Earth radius in cm

savefigures = False

files_directory = 'ACI/Model_Optimization/Clustered_P_R/All_params_random_weighted_targs150060_maxincl80/'
loadfiles_directory = 'ExoplanetsSysSim.jl-master/examples/clusters/' + files_directory
savefigures_directory = 'Clustering_Method_Figures/ExoplanetsSysSim/Power_law_r1_r2_sigma_r/Optimization_Plots/' + files_directory + 'Best_models/'

#files_directory = 'ACI/Model_Optimization/Non_clustered/Some9_params1_random_weighted_targs150060_maxincl80/'
#loadfiles_directory = 'ExoplanetsSysSim.jl-master/examples/clusters/' + files_directory
#savefigures_directory = 'Clustering_Method_Figures/ExoplanetsSysSim/Non_clustered/Optimization_Plots/' + files_directory + 'Best_models/'

subdirectory = 'Talk_Figures/'
model_name = 'ExoplanetsSysSim_Clustered_Model' #'ExoplanetsSysSim_Clustered_Model'





#To first read the number of simulated targets and bounds for the periods and radii:
with open(loadfiles_directory + 'periods1.out', 'r') as file:
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

#To read in the text file with all the best model parameters (just to get the first column of run numbers so that we can save the figures using the same numbering):
active_params_best_weighted_table = np.genfromtxt(loadfiles_directory + 'Active_params_best_weighted_all.txt', names=True)

res_ratios, res_width = [1.5, 2.0], 0.05 #NOTE: in the model, the near-resonant planets have period ratios between X and (1+w)*X where X = [2/1, 3/2, 4/3, 5/4] and w = 0.05!




##### To define functions which computes KS and AD distances:
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





##### To load and compute the exoplanet multiplicities, periods, and period ratios of the confirmed Kepler exoplanets:
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
        system_sort_i = np.argsort(system_P) #indices that would sort the periods of the planets in this system
        system_P = system_P[system_sort_i] #periods of all the planets in this system, sorted
        system_t_D = system_t_D[system_sort_i] #transit durations of all the planets in this system, sorted by period
        system_D = system_D[system_sort_i] #transit depths of all the planets in this system, sorted by period
        
        #To count the total number of planets in this system:
        M_confirmed.append(len(system_P))
        
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
P_confirmed = table_confirmed['P']
M_confirmed = np.array(M_confirmed)
R_confirmed = np.array(R_confirmed)
D_ratio_confirmed = np.array(D_ratio_confirmed)
xi_confirmed = np.array(xi_confirmed)
xi_res_confirmed = np.array(xi_res_confirmed)
xi_res32_confirmed = np.array(xi_res32_confirmed)
xi_res21_confirmed = np.array(xi_res21_confirmed)
xi_nonres_confirmed = np.array(xi_nonres_confirmed)





##### To load the files with the systems with observed planets and plot them:
'''
param_keys_all = [("num_targets_sim_pass_one", r'$N_{\rm stars,sim}$'),
                  ("max_incl_sys", r'$i_{\rm ref,max}$'),
                  ("log_rate_clusters", r'$\lambda_c$'),
                  ("max_clusters_in_sys", r'$N_{c,\rm max}$'),
                  ("power_law_P", r'$\alpha_P$'),
                  ("min_period", r'$P_{\rm min}$'),
                  ("max_period", r'$P_{\rm max}$'),
                  ("power_law_r1", r'$\alpha_{R1}$'),
                  ("power_law_r2", r'$\alpha_{R2}$'),
                  ("min_radius (R_earth)", r'$R_{p,\rm min}$ $(R_\oplus)$'),
                  ("max_radius (R_earth)", r'$R_{p,\rm max}$ $(R_\oplus)$'),
                  ("break_radius (R_earth)", r'$R_{p,\rm break}$ $(R_\oplus)$'),
                  ("sigma_incl", r'$\sigma_i$'),
                  ("sigma_incl_near_mmr", r'$\sigma_{i,\rm res}$'),
                  ("sigma_hk", r'$\sigma_e$'),
                  ("num_mutual_hill_radii", r'$\Delta_c$'),
                  ("mr_power_index", r'$\alpha_{mr}$'),
                  ("mr_max_mass (M_earth)", r'$M_{p,\rm max}$ $(M_\oplus)$')] #list of the symbols and names for all the model parameters; NOTE: although the params are named log rate of clusters and planets per cluster, we use the symbols and values for the rates
'''
#'''
param_keys_all = [("num_targets_sim_pass_one", r'$N_{\rm stars,sim}$'),
                  ("max_incl_sys", r'$i_{\rm ref,max}$'),
                  ("log_rate_clusters", r'$\lambda_c$'),
                  ("max_clusters_in_sys", r'$N_{c,\rm max}$'),
                  ("log_rate_planets_per_cluster", r'$\lambda_p$'),
                  ("max_planets_in_clusters", r'$N_{p,\rm max}$'),
                  ("power_law_P", r'$\alpha_P$'),
                  ("min_period", r'$P_{\rm min}$'),
                  ("max_period", r'$P_{\rm max}$'),
                  ("power_law_r1", r'$\alpha_{R1}$'),
                  ("power_law_r2", r'$\alpha_{R2}$'),
                  ("min_radius (R_earth)", r'$R_{p,\rm min}$ $(R_\oplus)$'),
                  ("max_radius (R_earth)", r'$R_{p,\rm max}$ $(R_\oplus)$'),
                  ("break_radius (R_earth)", r'$R_{p,\rm break}$ $(R_\oplus)$'),
                  ("sigma_incl", r'$\sigma_i$'),
                  ("sigma_incl_near_mmr", r'$\sigma_{i,\rm res}$'),
                  ("sigma_hk", r'$\sigma_e$'),
                  ("num_mutual_hill_radii", r'$\Delta_c$'),
                  ("mr_power_index", r'$\alpha_{mr}$'),
                  ("mr_max_mass (M_earth)", r'$M_{p,\rm max}$ $(M_\oplus)$'),
                  ("sigma_log_radius_in_cluster", r'$\sigma_R$'),
                  ("sigma_logperiod_per_pl_in_cluster", r'$\sigma_N$')] #list of the symbols and names for all the model parameters; NOTE: although the params are named log rate of clusters and planets per cluster, we use the symbols and values for the rates
#'''

#'''
for i in [1]: #active_params_best_weighted_table['run_number']
    run_number = int(i)

    #To read the simulation parameters from the file:
    param_vals_all = [] #list to be filled with the values of all the model parameters
    with open(loadfiles_directory + 'periods%s.out' % run_number, 'r') as file:
        for line in file:
            for i in range(len(param_keys_all)):
                chars = len(param_keys_all[i][0])
                if line[:3+chars] == '# ' + param_keys_all[i][0] + ':':
                    if param_keys_all[i][0][:3] == 'log':
                        param_vals_all.append(np.round(np.exp(float(line[4+chars:])), 4))
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
    radii_obs = np.sqrt(D_obs)*np.transpose([Rstar_obs])*(Rsun/Rearth) #array of planet radii, in Earth radii



    #To calculate the observed period ratios, period-normalized transit duration ratios, and transit depth ratios:
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





    ##### To compare the simulated observed distributions to the Kepler observed distributions using the K-S distance:

    #To compute the K-S distances and their positions, as well as additional statistics:
    delta_f = np.abs(len(P_obs)/float(N_sim) - len(P_confirmed)/float(N_Kep)) #absolute difference in the rates of observed planets per star

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



    #To plot the 'observed' distributions with the actual observed Kepler distributions:
    linewidth = 1
    if subdirectory == 'Talk_Figures/':
        linewidth = 3

    fig = plt.figure(figsize=(16,8))
    plot = GridSpec(4,2,left=0.075,bottom=0.1,right=0.95,top=0.95,wspace=0.15,hspace=0.4)

    #To print the parameter values:
    nrows = 8
    for i in range(len(param_keys_all)): #range(len(param_keys_all))
        plt.figtext(x=0.52+0.14*int(i/float(nrows)), y=0.95-0.025*(i%nrows), s=r'%s = %s' % (param_keys_all[i][1], param_vals_all[i]), fontsize=12)

    ax = plt.subplot(plot[0,0])
    x = Mtot_obs[Mtot_obs > 0]
    max_M = np.max((np.max(Mtot_obs), np.max(M_confirmed)))
    counts, bins = np.histogram(x, bins=max_M+1, range=(-0.5, max_M+0.5))
    bins_mid = (bins[:-1] + bins[1:])/2.
    plt.plot(bins_mid, counts/float(np.sum(counts)), 'o-', color='k', linewidth=linewidth, label='%s simulated systems' % len(x))
    counts, bins = np.histogram(M_confirmed, bins=bins)
    plt.plot(bins_mid, counts/float(np.sum(counts)), 'o--', color='k', alpha=0.2, label='%s Kepler systems' % len(M_confirmed))
    plt.gca().set_yscale("log")
    ax.tick_params(axis='both', labelsize=12)
    plt.xlim([1., max_M])
    if subdirectory == 'Talk_Figures/':
        plt.xlim([1., 8.])
    plt.xlabel(r'$M_{\rm tot}$', fontsize=12)
    plt.ylabel('Fraction', fontsize=12)
    plt.legend(loc='lower left', bbox_to_anchor=(0.01,0.01), frameon=False, ncol=1, fontsize=12) #show the legend
    plt.figtext(x=0.47, y=0.93, s=r'$\mathcal{D} = %1.4f$' % M_KS, ha='right', fontsize=12)

    ax = plt.subplot(plot[1,0])
    hist = plt.hist(P_obs_flat, bins=np.logspace(np.log10(P_min), np.log10(P_max), 101), histtype='step', weights=np.ones(len(P_obs_flat))/len(P_obs_flat), log=True, color='k', linewidth=linewidth, label=r'Simulated')
    plt.hist(P_confirmed, bins=hist[1], histtype='stepfilled', weights=np.ones(len(P_confirmed))/len(P_confirmed), log=True, color='k', alpha=0.2, label=r'Kepler')
    plt.gca().set_xscale("log")
    ax.tick_params(axis='both', labelsize=12)
    ax.set_xticks([3,10,30,100,300])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.xlim([P_min, 1.1*P_max])
    plt.ylim([10.**(-3.), 0.1])
    plt.xlabel(r'$P$ (days)', fontsize=12)
    plt.ylabel('Fraction', fontsize=12)
    plt.figtext(x=0.47, y=0.7, s=r'$\mathcal{D} = %1.4f$' % P_KS, ha='right', fontsize=12)

    ax = plt.subplot(plot[1,1])
    R_max_cut = 30. #upper cut-off for plotting period ratios; np.max(Rm_obs_flat)
    x = Rm_obs_flat[Rm_obs_flat < R_max_cut]
    hist = plt.hist(x, bins=np.logspace(np.log10(min(np.min(x),np.min(R_confirmed))), np.log10(R_max_cut), 101), histtype='step', weights=np.ones(len(x))/len(x), color='k', linewidth=linewidth, label='Simulated')
    x = R_confirmed[R_confirmed < R_max_cut]
    for i in range(len(res_ratios)):
        plt.axvline(x=res_ratios[i], linestyle=':', color='k')
    plt.hist(x, bins=hist[1], histtype='stepfilled', weights=np.ones(len(x))/len(x), color='k', alpha=0.2, label='Kepler')
    plt.gca().set_xscale("log")
    ax.tick_params(axis='both', labelsize=12)
    ax.set_xticks([1,2,3,4,5,10,20])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.xlim([1, R_max_cut])
    plt.xlabel(r'$P_{i+1}/P_i$', fontsize=12)
    plt.figtext(x=0.94, y=0.7, s=r'$\mathcal{D} = %1.4f$' % R_KS, ha='right', fontsize=12)

    ax = plt.subplot(plot[2,0])
    x = tdur_obs_flat
    hist = plt.hist(x, bins=np.linspace(min(np.min(tdur_obs_flat),np.min(t_D_confirmed*60.)), max(np.max(tdur_obs_flat),np.max(t_D_confirmed*60.)), 101), histtype='step', weights=np.ones(len(x))/len(x), color='k', linewidth=linewidth, label='Simulated')
    plt.hist(t_D_confirmed*60., bins=hist[1], histtype='stepfilled', weights=np.ones(len(t_D_confirmed))/len(t_D_confirmed), color='k', alpha=0.2, label=r'Kepler')
    ax.tick_params(axis='both', labelsize=12)
    plt.xlim([hist[1][0], hist[1][-1]])
    if subdirectory == 'Talk_Figures/':
        plt.xlim([0., 1500.])
        plt.ylim([0., 0.08])
    plt.xlabel(r'$t_{\rm dur}$ (mins)', fontsize=12)
    plt.ylabel('Fraction', fontsize=12)
    plt.figtext(x=0.47, y=0.47, s=r'$\mathcal{D} = %1.4f$' % tdur_KS, ha='right', fontsize=12)

    ax = plt.subplot(plot[2,1])
    x = np.log10(xi_obs_flat)
    hist = plt.hist(x, bins=np.linspace(-0.5, 0.5, 101), histtype='step', weights=np.ones(len(x))/len(x), color='k', linewidth=linewidth, label='All')
    plt.hist(np.log10(xi_confirmed), bins=hist[1], histtype='stepfilled', weights=np.ones(len(xi_confirmed))/float(len(xi_confirmed)), color='k', alpha=0.2)
    ax.tick_params(axis='both', labelsize=12)
    plt.xlim([hist[1][0], hist[1][-1]])
    if subdirectory == 'Talk_Figures/':
        plt.ylim([0., 0.1])
    plt.xlabel(r'$\log{\xi}$', fontsize=12)
    plt.figtext(x=0.94, y=0.47, s=r'$\mathcal{D} = %1.4f$' % logxi_KS, ha='right', fontsize=12)

    ax = plt.subplot(plot[3,0])
    x = D_obs_flat
    hist = plt.hist(x, bins=np.logspace(np.log10(min(np.min(x),np.min(D_confirmed))), np.log10(max(np.max(x),np.max(D_confirmed))), 101), histtype='step', weights=np.ones(len(x))/len(x), color='k', linewidth=linewidth, label='Simulated')
    plt.hist(D_confirmed, bins=hist[1], histtype='stepfilled', weights=np.ones(len(D_confirmed))/len(D_confirmed), color='k', alpha=0.2, label=r'Kepler')
    plt.gca().set_xscale("log")
    ax.tick_params(axis='both', labelsize=12)
    plt.xlim([hist[1][0], hist[1][-1]])
    if subdirectory == 'Talk_Figures/':
        plt.xlim([1e-5, 1e-2])
        plt.ylim([0., 0.04])
    plt.xlabel(r'$\delta$', fontsize=12)
    plt.ylabel('Fraction', fontsize=12)
    plt.figtext(x=0.47, y=0.24, s=r'$\mathcal{D} = %1.4f$' % D_KS, ha='right', fontsize=12)

    ax = plt.subplot(plot[3,1])
    x = D_ratio_obs_flat
    hist = plt.hist(x, bins=np.logspace(np.log10(min(np.min(x),np.min(D_ratio_confirmed))), np.log10(max(np.max(x),np.max(D_ratio_confirmed))), 101), histtype='step', weights=np.ones(len(x))/len(x), color='k', linewidth=linewidth, label='Simulated')
    plt.hist(D_ratio_confirmed, bins=hist[1], histtype='stepfilled', weights=np.ones(len(D_ratio_confirmed))/len(D_ratio_confirmed), color='k', alpha=0.2, label=r'Kepler')
    plt.gca().set_xscale("log")
    ax.tick_params(axis='both', labelsize=12)
    plt.xlim([hist[1][0], hist[1][-1]])
    if subdirectory == 'Talk_Figures/':
        plt.xlim([1e-2, 1e2])
        plt.ylim([0., 0.08])
    plt.xlabel(r'$\delta_{i+1}/\delta_i$', fontsize=12)
    plt.figtext(x=0.94, y=0.24, s=r'$\mathcal{D} = %1.4f$' % D_ratio_KS, ha='right', fontsize=12)

    if savefigures == True:
        plt.savefig(savefigures_directory + subdirectory + model_name + '_%s_observed_summary.pdf' % run_number)
    else:
        plt.show()
    plt.close()
#'''






