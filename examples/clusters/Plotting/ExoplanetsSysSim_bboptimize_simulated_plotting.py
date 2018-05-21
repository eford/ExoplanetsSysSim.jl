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





##### This module will be used to plot results of the optimization runs of our clustered model using bboptimize:

#To define some useful constants:
AU = 1.496*10.**13. #AU in cm
Msun = 1.989*10.**30. #Solar mass in kg
Rsun = 6.957*10.**10. #Solar radius in cm
Rearth = 6.371*10.**8. #Earth radius in cm

savefigures = False

loadfiles_directory = 'ACI/Model_Optimization/Clustered_P_R/Simulated_catalog_optimization/Some6_params_random_weighted_targs150060_maxincl80/' #All_params_random_weighted_targs150060_maxincl80/; Some10_params_random_weighted_targs150060_maxincl80/; Some6_params_random_weighted_targs150060_maxincl80/
weighted = True

savefigures_directory = 'Clustering_Method_Figures/ExoplanetsSysSim/Power_law_r1_r2_sigma_r/Optimization_Plots/' + loadfiles_directory

model_name = 'ExoplanetsSysSim_Clustered_Model_bboptimize' #'ExoplanetsSysSim_Clustered_Model_bboptimize'





##### To iterate through each of the optimization runs (files), and extract the results:

distances_names = [r'$|f_{\rm sim} - f_{\rm Kep}|$', 'Multiplicity', 'Period', 'Period ratio', 'Duration', 'xi', r'$\delta$', r'$\delta_{i+1}/\delta_i$'] #names of the distances; everything except the first entry is a KS distance of that distribution; [r'$|f_{\rm sim} - f_{\rm Kep}|$', 'Multiplicity', 'Period', 'Period ratio', 'Duration', 'xi', r'$\delta$', r'$\delta_{i+1}/\delta_i$']

mean_distances_perfect_all = [] #list to be filled with arrays of the distances of a good model compared to itself 20 times
std_distances_perfect_all = [] #list to be filled with arrays of the standard deviations of the distances of a good model compared to itself 20 times, which serves as the weights for the weighted distances
mean_weighted_distances_perfect_all = [] #list to be filled with arrays of the weighted distances (mean/std of distance) of a good model compared to itself 20 times
tolerance_all = [] #list to be filled with the tolerance (+/-) distance based on a good model compared to itself 20 times, for each run (this is either the std of the total distances or the total weighted distances)
time_perfect_all = [] #list to be filled with the elapsed times (s) for running a good model 20 times and comparing to itself, for each run

active_params_true_all = [] #list to be filled with arrays of the true values of the active parameters for each run (should be the same for all runs)
active_params_names_all = [] #list to be filled with arrays of the names of the active parameters for each run (should be the same for all runs)
active_params_bounds_all = [] #list to be filled with arrays of the search bounds of the active parameters for each run (should be the same for all runs)
active_params_start_all = [] #list to be filled with arrays of the starting values of the active parameters for each run
active_params_best_all = [] #list to be filled with arrays of the best values (lowest total distance) of the active parameters for each run
active_params_best_weighted_all = [] #list to be filled with arrays of the best values (lowest total weighted distance) of the active parameters for each run
distances_start_all = [] #list to be filled with arrays of the distances of the model with the starting active parameters compared to the Kepler sample, for each run
weighted_distances_start_all = [] #list to be filled with arrays of the weighted distances of the model with the starting active parameters compared to the Kepler sample, for each run
distances_best_all = [] #list to be filled with arrays of the distances of the model with the best values (i.e. summing to the lowest distance) compared to the Kepler sample, for each run
weighted_distances_best_all = [] #list to be filled with arrays of the weighted distances of the model with the best values (i.e. summing to the lowest weighted distance) compared to the Kepler sample, for each run
steps_best_all = [] #list to be filled with the number of model iterations to find the best active parameter values (lowest total distance) for each run
steps_best_weighted_all = [] #list to be filled with the number of model iterations to find the best active parameter values (lowest total weighted distance) for each run
steps_tot_all = [] #list to be filled with the number of total model iterations in the optimization procedure, for each run
time_optimization_all = [] #list to be filled with the elapsed times (s) for the full optimization procedure, for each run

active_params_steps_all = [] #list to be filled with arrays of the values of all the active parameters at every step (excluding starting values but including best values) for all the runs
distances_steps_all = [] #list to be filled with arrays of the distances of the model compared to the Kepler sample at all steps of the optimizations for all the runs
weighted_distances_steps_all = [] #list to be filled with arrays of the weighted distances of the model compared to the Kepler sample at all steps of the optimizations for all the runs

runs_started = 0
runs_finished = 0
for i in range(1,21): #range(1,21)
    with open(loadfiles_directory + 'Clustered_P_R_broken_R_simulated_optimization_random%s_targs150060_evals1000.txt' % i, 'r') as file:
        optim_lines = False #set to true once we start reading lines in the file that are the outputs of the optimization
        active_params_start = [] #will be replaced by the actual active parameter values if the file is not empty
        active_params_best = [] #will be replaced by the actual best active parameter values (lowest total distance) if the optimization progressed
        active_params_best_weighted = [] #will be replaced by the actual best active parameter values (lowest total weighted distance) if the optimization progressed
        distances_best = [1]*len(distances_names) #will be replaced with the best distances as the file is read, if the optimization progressed
        weighted_distances_best = [1e6]*len(distances_names) #will be replaced with the best weighted distances as the file is read, if the optimization progressed
        best_distance = sum(distances_best) #will be replaced with the best total distance if the optimization progressed
        best_fitness = sum(weighted_distances_best) #will be replaced with the best total weighted distance if the optimization progressed
        steps = 0 #will be a running count of the number of model iterations
        steps_best = steps #will be replaced by the number of the model iteration at which the best total distance was found
        steps_best_weighted = steps #will be replaced by the number of the model iteration at which the best total weighted distance was found
        for line in file:
            #For recording the preliminary runs of the model before optimizations:
            if line[0:13] == 'Active_params' and optim_lines == False and len(active_params_true_all) < i:
                active_params = [float(x) for x in line[16:-2].split(', ')]
                active_params_true_all.append(active_params)
            elif line[0:5] == 'Mean:':
                mean_distances_perfect_str, mean_distance_tot_perfect_str = line[7:-2].split('][')
                mean_distances_perfect, mean_distance_tot_perfect = [float(x) for x in mean_distances_perfect_str.split(', ')], float(mean_distance_tot_perfect_str)
                mean_distances_perfect_all.append(mean_distances_perfect)
            elif line[0:3] == 'Std':
                std_distances_perfect_str, std_distance_tot_perfect_str = line[7:-2].split('][')
                std_distances_perfect, std_distance_tot_perfect = [float(x) for x in std_distances_perfect_str.split(', ')], float(std_distance_tot_perfect_str)
                std_distances_perfect_all.append(std_distances_perfect)
            elif line[0:5] == 'Mean ':
                mean_weighted_distances_perfect_str, mean_weighted_distance_tot_perfect_str = line[22:-2].split('][')
                mean_weighted_distances_perfect, mean_weighted_distance_tot_perfect = [float(x) for x in mean_weighted_distances_perfect_str.split(', ')], float(mean_weighted_distance_tot_perfect_str)
                mean_weighted_distances_perfect_all.append(mean_weighted_distances_perfect)
            elif line[0:10] == '# Distance':
                if weighted == False:
                    target_str, tolerance_str = line[57:-1].split(' +/- ')
                    tolerance = float(tolerance_str)
                    tolerance_all.append(tolerance)
            elif line[0:10] == '# Weighted':
                if weighted == True:
                    target_str, tolerance_str = line[66:-1].split(' +/- ')
                    tolerance = float(tolerance_str)
                    tolerance_all.append(tolerance)
            elif line[0:9] == '# elapsed' and optim_lines == False:
                time_perfect_all.append(float(line[16:-8]))
            elif line[0:19] == '# Active parameters':
                active_params_names = line[29:-3].split('", "')
                active_params_names_all.append(active_params_names)
            
            #For recording the results of the optimizations:
            elif line[0:7] == '# Start':
                active_params_start = [float(x) for x in line[37:-2].split(', ')]
            elif line[0:7] == '# Optim':
                runs_started += 1
                active_params_bounds = [(float(x.split(', ')[0]), float(x.split(', ')[1])) for x in line[73:-3].split('), (')]
                optim_lines = True
            elif line[0:13] == 'Active_params' and optim_lines == True:
                active_params = [float(x) for x in line[16:-2].split(', ')]
                active_params_steps_all.append(active_params)
            elif line[0:5] == 'Dist:' and optim_lines == True:
                steps += 1
                distances_str, distance_tot_str = line[7:-2].split('][')
                distances, distance_tot = [float(x) for x in distances_str.split(', ')], float(distance_tot_str)
                distances_steps_all.append(distances)
                if distance_tot < sum(distances_best):
                    distances_best = distances
                    best_distance = sum(distances_best)
                    if weighted == False:
                        best_fitness = sum(distances_best) #If optimizer used sum of distances
                    active_params_best = active_params
                    steps_best = steps
                if steps == 1:
                    distances_start_all.append(distances)
            elif line[0:13] == 'Dist_weighted' and optim_lines == True:
                weighted_distances_str, weighted_distance_tot_str = line[16:-2].split('][')
                weighted_distances, weighted_distance_tot = [float(x) for x in weighted_distances_str.split(', ')], float(weighted_distance_tot_str)
                weighted_distances_steps_all.append(weighted_distances)
                if weighted_distance_tot < sum(weighted_distances_best):
                    weighted_distances_best = weighted_distances
                    if weighted == True:
                        best_fitness = sum(weighted_distances_best) #If optimizer used weighted distances
                    active_params_best_weighted = active_params
                    steps_best_weighted = steps
                if steps == 1:
                    weighted_distances_start_all.append(weighted_distances)
            elif line[0:14] == '# best_fitness':
                runs_finished += 1
                best_fitness_end = float(line[16:-2])
                #print i
                #print i, best_fitness_end, sum(weighted_distances_best), best_fitness-sum(weighted_distances_best)
            elif line[0:9] == '# elapsed' and optim_lines == True:
                time_optimization_all.append(float(line[16:-8]))
    
        print i, optim_lines, best_fitness, len(active_params_steps_all), len(distances_steps_all), len(weighted_distances_steps_all)
        if best_fitness < 8.*1e6 and optim_lines == True: #only keep the runs in which the optimization actually progressed (i.e. to discard killed or otherwise faulty runs)
            active_params_bounds_all.append(active_params_bounds)
            active_params_start_all.append(active_params_start)
            active_params_best_all.append(active_params_best)
            active_params_best_weighted_all.append(active_params_best_weighted)
            distances_best_all.append(distances_best)
            weighted_distances_best_all.append(weighted_distances_best)
            steps_best_all.append(steps_best)
            steps_best_weighted_all.append(steps_best_weighted)
            steps_tot_all.append(steps)

print 'Runs successfully started (and not killed): ', runs_started #runs killed because of the wall time are not counted here because they have their output files emptied
print 'Runs successfully finished (reached max iterations or target fitness): ', runs_finished #runs not counted here are ones killed either because of the wall time, or because of bus error

mean_distances_perfect_all = np.array(mean_distances_perfect_all)
mean_distance_tot_perfect_all = np.sum(mean_distances_perfect_all, axis=1)
std_distances_perfect_all = np.array(std_distances_perfect_all)
mean_weighted_distances_perfect_all = np.array(mean_weighted_distances_perfect_all)
mean_weighted_distance_tot_perfect_all = np.sum(mean_weighted_distances_perfect_all, axis=1)
tolerance_all = np.array(tolerance_all)
time_perfect_all = np.array(time_perfect_all)

active_params_true_all = np.array(active_params_true_all)
active_params_names_all = np.array(active_params_names_all)
active_params_bounds_all = np.array(active_params_bounds_all)
active_params_start_all = np.array(active_params_start_all)
active_params_best_all = np.array(active_params_best_all)
active_params_best_weighted_all = np.array(active_params_best_weighted_all)
distances_start_all = np.array(distances_start_all)
weighted_distances_start_all = np.array(weighted_distances_start_all)
distances_best_all = np.array(distances_best_all)
weighted_distances_best_all = np.array(weighted_distances_best_all)

distance_tot_start_all = np.sum(distances_start_all, axis=1)
if weighted == True:
    weighted_distance_tot_start_all = np.sum(weighted_distances_start_all, axis=1)
distance_tot_best_all = np.sum(distances_best_all, axis=1)
if weighted == True:
    weighted_distance_tot_best_all = np.sum(weighted_distances_best_all, axis=1)

steps_best_all = np.array(steps_best_all)
steps_best_weighted_all = np.array(steps_best_weighted_all)
steps_tot_all = np.array(steps_tot_all)
time_optimization_all = np.array(time_optimization_all)

active_params_steps_all = np.array(active_params_steps_all)
distances_steps_all = np.array(distances_steps_all)
weighted_distances_steps_all = np.array(weighted_distances_steps_all)
distance_tot_steps_all = np.sum(distances_steps_all, axis=1)
if weighted == True:
    weighted_distance_tot_steps_all = np.sum(weighted_distances_steps_all, axis=1)

if weighted == False:
    convergence_all = distance_tot_best_all < (mean_distance_tot_perfect_all + tolerance_all) #array of booleans specifying whether each run converged to better than the target fitness +/- tolerance or not
elif weighted == True:
    convergence_all = weighted_distance_tot_best_all < (mean_weighted_distance_tot_perfect_all + tolerance_all) #array of booleans specifying whether each run converged to better than the target fitness +/- tolerance or not





##### To make 1D plots of the distances vs. each parameter:

#'''
N_steps_sample = min(1000, len(active_params_steps_all))
i_steps_sample = np.random.choice(np.arange(len(active_params_steps_all)), N_steps_sample, replace=False) #array of indices of a sample of optimization steps to be plotted

for i,param in enumerate(active_params_names):
    fig = plt.figure(figsize=(16,8))
    plot = GridSpec(len(distances_names)+1,1,left=0.075,bottom=0.115,right=0.875,top=0.925,wspace=0.25,hspace=0)
    ax = plt.subplot(plot[0,0])
    i_sortx = np.argsort(active_params_best_all[:,i]) #array of indices that would sort the array based on the current active parameter
    plt.axvline(x=active_params_true_all[0][i], linewidth=5, color='g') #to plot the true value
    plt.plot(active_params_best_all[i_sortx,i], distance_tot_best_all[i_sortx], 'o-', color='r', label='Total') #to plot the best values for each run
    plt.scatter(active_params_steps_all[i_steps_sample,i], distance_tot_steps_all[i_steps_sample], marker='.', color='k', alpha=0.1) #to plot the values for a sample of all the steps of all the runs
    plt.xlim(active_params_bounds_all[0][i])
    ax.set_xticks([])
    plt.yticks([np.round(1.1*np.min(distance_tot_steps_all),4), np.round(0.9*np.max(distance_tot_steps_all),4)])
    ax.tick_params(axis='both', labelsize=12)
    plt.legend(loc='center left', bbox_to_anchor=(1.,0.5), ncol=1, fontsize=12)
    for j,dist in enumerate(distances_names):
        ax = plt.subplot(plot[j+1,0])
        i_sortx = np.argsort(active_params_best_all[:,i]) #array of indices that would sort the array based on the current active parameter
        plt.axvline(x=active_params_true_all[0][i], linewidth=5, color='g') #to plot the true value
        plt.plot(active_params_best_all[i_sortx,i], distances_best_all[i_sortx,j], 'o-', color='r', label=dist)
        plt.scatter(active_params_steps_all[i_steps_sample,i], distances_steps_all[i_steps_sample,j], marker='.', color='k', alpha=0.1)
        plt.xlim(active_params_bounds_all[0][i])
        if j != len(distances_names)-1:
            ax.set_xticks([])
        plt.yticks([np.round(1.1*np.min(distances_steps_all[:,j]),4), np.round(0.9*np.max(distances_steps_all[:,j]),4)])
        ax.tick_params(axis='both', labelsize=12)
        if j == 3:
            plt.ylabel('Distance (best)', fontsize=20)
        plt.legend(loc='center left', bbox_to_anchor=(1.,0.5), ncol=1, fontsize=12)
    plt.xlabel(param + ' (best)', fontsize=20)

    if savefigures == True:
        plt.savefig(savefigures_directory + model_name + '_' + param + '.pdf')
    else:
        plt.show()
    plt.close()
#'''




##### To make 2D plots of various pairs of parameters:

#'''
#active_params_pairs = [("log_rate_clusters", "log_rate_planets_per_cluster"), ("power_law_r1", "power_law_r2"), ("sigma_hk", "sigma_incl"), ("num_mutual_hill_radii", "sigma_logperiod_per_pl_in_cluster"), ("break_radius", "sigma_log_radius_in_cluster"), ("break_radius", "power_law_r2"), ("mr_power_index", "power_law_P"), ("sigma_incl", "sigma_incl_near_mmr")] #for all (13) active parameters (clustered model)
#active_params_pairs = [("log_rate_clusters", "log_rate_planets_per_cluster"), ("power_law_r1", "power_law_r2"), ("sigma_hk", "sigma_incl"), ("power_law_P", "sigma_logperiod_per_pl_in_cluster"), ("break_radius", "sigma_log_radius_in_cluster"), ("break_radius", "power_law_r2")] #for some (10) active parameters (clustered model)
active_params_pairs = [("log_rate_clusters", "log_rate_planets_per_cluster"), ("power_law_P", "power_law_r1"), ("power_law_r1", "power_law_r2"), ("break_radius", "power_law_r1"), ("break_radius", "power_law_r2")] #for some (6) active parameters (clustered model)

for i,pair in enumerate(active_params_pairs):
    i_x, i_y = np.where(np.array(active_params_names) == pair[1])[0][0], np.where(np.array(active_params_names) == pair[0])[0][0]

    #To plot the total distances:
    fig = plt.figure(figsize=(16,8))
    plot = GridSpec(1,1,left=0.1,bottom=0.115,right=0.95,top=0.925,wspace=0.25,hspace=0)
    ax = plt.subplot(plot[0,0])
    plt.scatter(active_params_true_all[0][i_x], active_params_true_all[0][i_y], marker='*', c='k', s=500, alpha=1) #true values
    best_scatter = plt.scatter(active_params_best_all[:,i_x], active_params_best_all[:,i_y], marker='.', c=distance_tot_best_all, s=500, alpha=1) #best values for each run
    plt.scatter(active_params_start_all[:,i_x], active_params_start_all[:,i_y], marker='.', c=distance_tot_start_all, s=100, alpha=1) #starting values for each run ### facecolors='none', edgecolors='k'
    #plt.scatter(active_params_steps_all[:,i_x], active_params_steps_all[:,i_y], marker='.', c=distance_tot_steps_all, s=50, alpha=0.5) #all values at each step of each run
    for j in range(len(active_params_start_all)):
        plt.plot([active_params_start_all[j,i_x],active_params_best_all[j,i_x]], [active_params_start_all[j,i_y],active_params_best_all[j,i_y]], '--', color='k') #to plot a line connecting the starting values to the best values
    #plt.axis('equal')
    plt.xlim(active_params_bounds_all[0][i_x])
    plt.ylim(active_params_bounds_all[0][i_y])
    ax.tick_params(axis='both', labelsize=20)
    plt.xlabel(active_params_names[i_x], fontsize=20)
    plt.ylabel(active_params_names[i_y], fontsize=20)
    plt.colorbar(best_scatter)
    #plt.legend(loc='center left', bbox_to_anchor=(1.,0.5), ncol=1, fontsize=12)

    if savefigures == True:
        plt.savefig(savefigures_directory + model_name + '_' + pair[0] + '_' + pair[1] + '_dist.pdf')
    else:
        plt.show()
    plt.close()

    #To plot the individual distances:
    fig = plt.figure(figsize=(16,8))
    plot = GridSpec(2,4,left=0.075,bottom=0.115,right=0.95,top=0.925,wspace=0.25,hspace=0.25)
    plot_rows = [0,0,0,0,1,1,1,1]
    plot_cols = [0,1,2,3,0,1,2,3]
    for j,dist in enumerate(distances_names): #for the individual distances
        ax = plt.subplot(plot[plot_rows[j],plot_cols[j]])
        plt.title(dist, fontsize=12)
        plt.scatter(active_params_true_all[0][i_x], active_params_true_all[0][i_y], marker='*', c='k', s=200, alpha=1) #true values
        best_scatter = plt.scatter(active_params_best_all[:,i_x], active_params_best_all[:,i_y], marker='.', c=distances_best_all[:,j], s=200, alpha=1) #best values for each run
        #plt.scatter(active_params_steps_all[:,i_x], active_params_steps_all[:,i_y], marker='.', c=distances_steps_all[:,j], s=50, alpha=0.5) #all values at each step of each run
        plt.xlim(active_params_bounds_all[0][i_x])
        plt.ylim(active_params_bounds_all[0][i_y])
        ax.tick_params(axis='both', labelsize=12)
        plt.colorbar(best_scatter)
    fig.text(0.5, 0.05, active_params_names[i_x], ha='center', fontsize=20)
    fig.text(0.025, 0.5, active_params_names[i_y], va='center', rotation='vertical', fontsize=20)
    #plt.legend(loc='center left', bbox_to_anchor=(1.,0.5), ncol=1, fontsize=12)

    if savefigures == True:
        plt.savefig(savefigures_directory + model_name + '_' + pair[0] + '_' + pair[1] + '_dists.pdf')
    else:
        plt.show()
    plt.close()
#'''

