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
AU = 1.496*10.**13. #AU in cm
Msun = 1.989*10.**30. #Solar mass in kg
Rsun = 6.957*10.**10. #Solar radius in cm
Rearth = 6.371*10.**8. #Earth radius in cm

#CHECK THESE NUMBERS!!!
Nstars_Kepler = 150061.
Nstars_sim = 150061.
Nplanets_Kepler = 2712.
cos_factor = np.cos(80.*np.pi/180.) #####

savefigures = False
savefigures_directory = 'Clustering_Method_Figures/ExoplanetsSysSim/'
model_name = 'ExoplanetsSysSim_Clustered_Model'





##### To load arrays/grids with model fit statistics (including KS distances):

#NOTE: in Julia, saving a 2D array into a file actually transposes the array! Thus the naming '...A_B_grid' actually refers to the varied parameters on the x (cols) and y (rows) of the saved arrays respectively; we will transpose them back after loading so that A refers to the parameter along y (rows) and B refers to the parameter along x (cols)

'''
Model_stats_P_r_grid = [] #list to be filled with grids (2D arrays) of model statistics
with open('Clustering_Method_Figures/ExoplanetsSysSim/Power_law_r/Model_stats_P_r_grid.out', 'r') as file: #ExoplanetsSysSim.jl-master/examples/clusters/
    for line in file:
        if line[0] != '#':
            grid = [] #list to be filled with lists (rows of the grid)
            line = line[1:-2]
            line_per_row = line.split('; ')
            #print len(line_per_sys)
            for x in line_per_row:
                run_list = x.split()
                run_list = [float(i) for i in run_list]
                grid.append(run_list)
            grid = np.transpose(np.array(grid))
            Model_stats_P_r_grid.append(grid)
            print np.shape(grid)
Model_stats_P_r_grid = np.array(Model_stats_P_r_grid)
'''

#'''
Model_stats_Nhill_sigma_grid = [] #list to be filled with grids (2D arrays) of model statistics
with open('Clustering_Method_Figures/ExoplanetsSysSim/Power_law_r1_r2_sigma_r/Model_stats_Nhill_sigmaP_grid.out', 'r') as file:
    for line in file:
        if line[0] != '#':
            grid = [] #list to be filled with lists (rows of the grid)
            line = line[1:-2]
            line_per_row = line.split('; ')
            #print len(line_per_sys)
            for x in line_per_row:
                run_list = x.split()
                run_list = [float(i) for i in run_list]
                grid.append(run_list)
            grid = np.transpose(np.array(grid))
            Model_stats_Nhill_sigma_grid.append(grid)
            print np.shape(grid)
Model_stats_Nhill_sigma_grid = np.array(Model_stats_Nhill_sigma_grid)
#'''
    
#'''
Model_stats_sigma_incl_sigma_hk_grid = [] #list to be filled with grids (2D arrays) of model statistics
with open('Clustering_Method_Figures/ExoplanetsSysSim/Power_law_r1_r2_sigma_r/Model_stats_sigma_incl_sigma_hk_grid.out', 'r') as file:
    for line in file:
        if line[0] != '#':
            grid = [] #list to be filled with lists (rows of the grid)
            line = line[1:-2]
            line_per_row = line.split('; ')
            #print len(line_per_sys)
            for x in line_per_row:
                run_list = x.split()
                run_list = [float(i) for i in run_list]
                grid.append(run_list)
            grid = np.transpose(np.array(grid))
            Model_stats_sigma_incl_sigma_hk_grid.append(grid)
            print np.shape(grid)
Model_stats_sigma_incl_sigma_hk_grid = np.array(Model_stats_sigma_incl_sigma_hk_grid)
#'''

'''
Model_stats_max_radius_mr_grid = [] #list to be filled with grids (2D arrays) of model statistics
with open('Clustering_Method_Figures/ExoplanetsSysSim/Power_law_r/Model_stats_max_radius_mr_grid.out', 'r') as file:
    for line in file:
        if line[0] != '#':
            grid = [] #list to be filled with lists (rows of the grid)
            line = line[1:-2]
            line_per_row = line.split('; ')
            #print len(line_per_sys)
            for x in line_per_row:
                run_list = x.split()
                run_list = [float(i) for i in run_list]
                grid.append(run_list)
            grid = np.transpose(np.array(grid))
            Model_stats_max_radius_mr_grid.append(grid)
            print np.shape(grid)
Model_stats_max_radius_mr_grid = np.array(Model_stats_max_radius_mr_grid)
'''

#'''
Model_stats_mr_sigmaR_grid = [] #list to be filled with grids (2D arrays) of model statistics
with open('Clustering_Method_Figures/ExoplanetsSysSim/Power_law_r1_r2_sigma_r/Model_stats_mr_sigmaR_grid.out', 'r') as file:
    for line in file:
        if line[0] != '#':
            grid = [] #list to be filled with lists (rows of the grid)
            line = line[1:-2]
            line_per_row = line.split('; ')
            #print len(line_per_sys)
            for x in line_per_row:
                run_list = x.split()
                run_list = [float(i) for i in run_list]
                grid.append(run_list)
            grid = np.transpose(np.array(grid))
            Model_stats_mr_sigmaR_grid.append(grid)
            print np.shape(grid)
Model_stats_mr_sigmaR_grid = np.array(Model_stats_mr_sigmaR_grid)
#'''
    
#'''
Model_stats_Nc_Np_grid = [] #list to be filled with grids (2D arrays) of model statistics
with open('Clustering_Method_Figures/ExoplanetsSysSim/Power_law_r1_r2_sigma_r/Model_stats_Nc_Np_grid.out', 'r') as file:
    for line in file:
        if line[0] != '#':
            grid = [] #list to be filled with lists (rows of the grid)
            line = line[1:-2]
            line_per_row = line.split('; ')
            #print len(line_per_sys)
            for x in line_per_row:
                run_list = x.split()
                run_list = [float(i) for i in run_list]
                grid.append(run_list)
            grid = np.transpose(np.array(grid))
            Model_stats_Nc_Np_grid.append(grid)
            print np.shape(grid)
Model_stats_Nc_Np_grid = np.array(Model_stats_Nc_Np_grid)
#'''

#'''
Model_stats_r1_r2_grid = [] #list to be filled with grids (2D arrays) of model statistics
with open('Clustering_Method_Figures/ExoplanetsSysSim/Power_law_r1_r2_sigma_r/Model_stats_r1_r2_rbreak2_grid.out', 'r') as file:
    for line in file:
        if line[0] != '#':
            grid = [] #list to be filled with lists (rows of the grid)
            line = line[1:-2]
            line_per_row = line.split('; ')
            #print len(line_per_sys)
            for x in line_per_row:
                run_list = x.split()
                run_list = [float(i) for i in run_list]
                grid.append(run_list)
            grid = np.transpose(np.array(grid))
            Model_stats_r1_r2_grid.append(grid)
            print np.shape(grid)
Model_stats_r1_r2_grid = np.array(Model_stats_r1_r2_grid)
#'''

'''
Model_stats_r1_r2_rbreak3_grid = [] #list to be filled with grids (2D arrays) of model statistics
with open('Clustering_Method_Figures/ExoplanetsSysSim/Power_law_r1_r2/Model_stats_r1_r2_rbreak3_grid.out', 'r') as file:
    for line in file:
        if line[0] != '#':
            grid = [] #list to be filled with lists (rows of the grid)
            line = line[1:-2]
            line_per_row = line.split('; ')
            #print len(line_per_sys)
            for x in line_per_row:
                run_list = x.split()
                run_list = [float(i) for i in run_list]
                grid.append(run_list)
            grid = np.transpose(np.array(grid))
            Model_stats_r1_r2_rbreak3_grid.append(grid)
            print np.shape(grid)
Model_stats_r1_r2_rbreak3_grid = np.array(Model_stats_r1_r2_rbreak3_grid)
'''

'''
Model_stats_r_sigma_r_grid = [] #list to be filled with grids (2D arrays) of model statistics
with open('Clustering_Method_Figures/ExoplanetsSysSim/Power_law_r_sigma_r/Model_stats_r_sigma_r_grid.out', 'r') as file:
    for line in file:
        if line[0] != '#':
            grid = [] #list to be filled with lists (rows of the grid)
            line = line[1:-2]
            line_per_row = line.split('; ')
            #print len(line_per_sys)
            for x in line_per_row:
                run_list = x.split()
                run_list = [float(i) for i in run_list]
                grid.append(run_list)
            grid = np.transpose(np.array(grid))
            Model_stats_r_sigma_r_grid.append(grid)
            print np.shape(grid)
Model_stats_r_sigma_r_grid = np.array(Model_stats_r_sigma_r_grid)
'''



Model_stats_fields = [r'$f_{\rm sim}/f_{\rm Kepler}$', 'K-S multiplicities', 'K-S periods', 'K-S period ratios', 'K-S transit durations', 'K-S xi', 'K-S transit depths', 'K-S transit depth ratios'] #names of the fields associated with the 'Model_stats_...' arrays

'''
##### To plot the 2D arrays of model fit statistics (N_obs, KS multiplicities, KS periods, KS period ratios, KS transit depths, KS xi) on the parameter space grids:

#To plot the 2D arrays of model fit statistics on the grid of 'power_law_P' vs. 'power_law_r':
power_law_P_range = np.linspace(-0.5, 0.5, 11)
power_law_r_range = np.linspace(-3., -1.2, 11)

d_y = np.diff(power_law_P_range)[0]
d_x = np.diff(power_law_r_range)[0]
xmin, xmax, ymin, ymax = power_law_r_range[0] - d_x/2., power_law_r_range[-1] + d_x/2., power_law_P_range[0] - d_y/2., power_law_P_range[-1] + d_y/2.
xy_ratio = (xmax - xmin)/(ymax - ymin)

fig = plt.figure(figsize=(16,8))
plot = GridSpec(2,3,left=0.075,bottom=0.115,right=0.95,top=0.925,wspace=0.15,hspace=0.4)

plot_rows = [0,0,0,1,1,1]
plot_cols = [0,1,2,0,1,2]
for i,name in enumerate(Model_stats_fields):
    ax = plt.subplot(plot[plot_rows[i],plot_cols[i]])
    plt.title(name, fontsize=20)
    data_array = Model_stats_P_r_grid[i]
    if i==0:
        data_array = (data_array/Nstars_sim)/(Nplanets_Kepler/Nstars_Kepler)
    plt.imshow(data_array, aspect=xy_ratio, origin='lower', extent=(xmin, xmax, ymin, ymax))
    plt.colorbar()
    plt.xlabel(r'power_law_r', fontsize=20)
    plt.ylabel(r'power_law_P', fontsize=20)

if savefigures == True:
    plt.savefig(savefigures_directory + 'Power_law_r/' + model_name + '_P_r_stats.pdf')
else:
    plt.show()
plt.close()
'''

'''
#To plot the 2D arrays of model fit statistics on the grid of 'num_mutual_hill_radii' vs. 'sigma_logperiod_per_pl_in_cluster':
num_mutual_hill_radii_range = np.linspace(5., 20., 11)
sigma_logperiod_per_pl_in_cluster_range = np.linspace(0.05, 0.5, 11)

d_y = np.diff(num_mutual_hill_radii_range)[0]
d_x = np.diff(sigma_logperiod_per_pl_in_cluster_range)[0]
xmin, xmax, ymin, ymax = sigma_logperiod_per_pl_in_cluster_range[0] - d_x/2., sigma_logperiod_per_pl_in_cluster_range[-1] + d_x/2., num_mutual_hill_radii_range[0] - d_y/2., num_mutual_hill_radii_range[-1] + d_y/2.
xy_ratio = (xmax - xmin)/(ymax - ymin)

fig = plt.figure(figsize=(16,8))
plot = GridSpec(2,4,left=0.075,bottom=0.115,right=0.95,top=0.925,wspace=0.15,hspace=0.4)

plot_rows = [0,0,0,0,1,1,1,1]
plot_cols = [0,1,2,3,0,1,2,3]
for i,name in enumerate(Model_stats_fields):
    ax = plt.subplot(plot[plot_rows[i],plot_cols[i]])
    plt.title(name, fontsize=20)
    data_array = Model_stats_Nhill_sigma_grid[i]
    if i==0:
        data_array = (data_array/Nstars_sim)/(Nplanets_Kepler/Nstars_Kepler)
    plt.imshow(data_array, aspect=xy_ratio, origin='lower', extent=(xmin, xmax, ymin, ymax))
    plt.colorbar()
    plt.xlabel(r'sigma_logperiod_per_pl_in_cluster', fontsize=20)
    plt.ylabel(r'num_mutual_hill_radii', fontsize=20)

if savefigures == True:
    plt.savefig(savefigures_directory + 'Power_law_r1_r2_sigma_r/' + model_name + '_Nhill_sigma_stats.pdf')
else:
    plt.show()
plt.close()
'''

'''
#To plot the 2D arrays of model fit statistics on the grid of 'sigma_incl' vs. 'sigma_hk':
sigma_incl_range = np.linspace(0., 5., 11)
sigma_hk_range = np.linspace(0., 0.2, 11)

d_y = np.diff(sigma_incl_range)[0]
d_x = np.diff(sigma_hk_range)[0]
xmin, xmax, ymin, ymax = sigma_hk_range[0] - d_x/2., sigma_hk_range[-1] + d_x/2., sigma_incl_range[0] - d_y/2., sigma_incl_range[-1] + d_y/2.
xy_ratio = (xmax - xmin)/(ymax - ymin)

fig = plt.figure(figsize=(16,8))
plot = GridSpec(2,4,left=0.075,bottom=0.115,right=0.95,top=0.925,wspace=0.15,hspace=0.4)

plot_rows = [0,0,0,0,1,1,1,1]
plot_cols = [0,1,2,3,0,1,2,3]
for i,name in enumerate(Model_stats_fields):
    ax = plt.subplot(plot[plot_rows[i],plot_cols[i]])
    plt.title(name, fontsize=20)
    data_array = Model_stats_sigma_incl_sigma_hk_grid[i]
    if i==0:
        data_array = (data_array/Nstars_sim)/(Nplanets_Kepler/Nstars_Kepler)
    plt.imshow(data_array, aspect=xy_ratio, origin='lower', extent=(xmin, xmax, ymin, ymax))
    plt.colorbar()
    plt.xlabel(r'sigma_hk', fontsize=20)
    plt.ylabel(r'sigma_incl', fontsize=20)

if savefigures == True:
    plt.savefig(savefigures_directory + 'Power_law_r1_r2_sigma_r/' + model_name + '_sigma_incl_sigma_hk_stats.pdf')
else:
    plt.show()
plt.close()
'''

'''
#To plot the 2D arrays of model fit statistics on the grid of 'max_radius' vs. 'mr_power_law':
max_radius_range = np.linspace(2., 12., 11)
mr_range = np.linspace(1.5, 3.5, 11)

d_y = np.diff(max_radius_range)[0]
d_x = np.diff(mr_range)[0]
xmin, xmax, ymin, ymax = mr_range[0] - d_x/2., mr_range[-1] + d_x/2., max_radius_range[0] - d_y/2., max_radius_range[-1] + d_y/2.
xy_ratio = (xmax - xmin)/(ymax - ymin)

fig = plt.figure(figsize=(16,8))
plot = GridSpec(2,3,left=0.075,bottom=0.115,right=0.95,top=0.925,wspace=0.15,hspace=0.4)

plot_rows = [0,0,0,1,1,1]
plot_cols = [0,1,2,0,1,2]
for i,name in enumerate(Model_stats_fields):
    ax = plt.subplot(plot[plot_rows[i],plot_cols[i]])
    plt.title(name, fontsize=20)
    data_array = Model_stats_max_radius_mr_grid[i]
    if i==0:
        data_array = (data_array/Nstars_sim)/(Nplanets_Kepler/Nstars_Kepler)
    plt.imshow(data_array, aspect=xy_ratio, origin='lower', extent=(xmin, xmax, ymin, ymax))
    plt.colorbar()
    plt.xlabel(r'mr_power_law', fontsize=20)
    plt.ylabel(r'max_radius (R_earth)', fontsize=20)

if savefigures == True:
    plt.savefig(savefigures_directory + 'Power_law_r/' + model_name + '_max_radius_mr_stats.pdf')
else:
    plt.show()
plt.close()
'''

'''
#To plot the 2D arrays of model fit statistics on the grid of 'mr_power_law' vs. 'sigma_log_radius_in_cluster':
mr_range = np.linspace(1.5, 3.5, 11)
sigmaR_range = np.linspace(0.1, 1.0, 11)

d_y = np.diff(mr_range)[0]
d_x = np.diff(sigmaR_range)[0]
xmin, xmax, ymin, ymax = sigmaR_range[0] - d_x/2., sigmaR_range[-1] + d_x/2., mr_range[0] - d_y/2., mr_range[-1] + d_y/2.
xy_ratio = (xmax - xmin)/(ymax - ymin)

fig = plt.figure(figsize=(16,8))
plot = GridSpec(2,4,left=0.075,bottom=0.115,right=0.95,top=0.925,wspace=0.15,hspace=0.4)

plot_rows = [0,0,0,0,1,1,1,1]
plot_cols = [0,1,2,3,0,1,2,3]
for i,name in enumerate(Model_stats_fields):
    ax = plt.subplot(plot[plot_rows[i],plot_cols[i]])
    plt.title(name, fontsize=20)
    data_array = Model_stats_mr_sigmaR_grid[i]
    if i==0:
        data_array = (data_array/Nstars_sim)/(Nplanets_Kepler/Nstars_Kepler)
    plt.imshow(data_array, aspect=xy_ratio, origin='lower', extent=(xmin, xmax, ymin, ymax))
    plt.colorbar()
    plt.xlabel(r'sigma_log_radius_in_cluster', fontsize=20)
    plt.ylabel(r'mr_power_law', fontsize=20)

if savefigures == True:
    plt.savefig(savefigures_directory + 'Power_law_r1_r2_sigma_r/' + model_name + '_mr_sigmaR_stats.pdf')
else:
    plt.show()
plt.close()
'''

'''
#To plot the 2D arrays of model fit statistics on the grid of 'log_rate_clusters' vs. 'log_rate_planets_per_cluster':
Model_stats_fields = [r'$f_{\rm sim}/f_{\rm Kepler}$', 'Number of observed singles', 'K-S multiplicities', 'K-S periods', 'K-S period ratios', 'K-S transit durations', 'K-S xi', 'K-S transit depths', 'K-S transit depth ratios'] #names of the fields associated with the 'Model_stats_...' arrays

Nc_range = np.linspace(1., 5., 11)
Np_range = np.linspace(1., 5., 11)

d_y = np.diff(Nc_range)[0]
d_x = np.diff(Np_range)[0]
xmin, xmax, ymin, ymax = Np_range[0] - d_x/2., Np_range[-1] + d_x/2., Nc_range[0] - d_y/2., Nc_range[-1] + d_y/2.
xy_ratio = (xmax - xmin)/(ymax - ymin)

fig = plt.figure(figsize=(16,8))
plot = GridSpec(3,3,left=0.075,bottom=0.115,right=0.95,top=0.925,wspace=0.25,hspace=0.4)

plot_rows = [0,0,0,1,1,1,2,2,2]
plot_cols = [0,1,2,0,1,2,0,1,2]
for i,name in enumerate(Model_stats_fields):
    ax = plt.subplot(plot[plot_rows[i],plot_cols[i]])
    plt.title(name, fontsize=20)
    data_array = Model_stats_Nc_Np_grid[i]
    if i==0:
        data_array = (data_array/Nstars_sim)/(Nplanets_Kepler/Nstars_Kepler)
    plt.imshow(data_array, aspect=xy_ratio, origin='lower', extent=(xmin, xmax, ymin, ymax))
    plt.colorbar()
    plt.xlabel(r'rate_planets_per_cluster', fontsize=20)
    plt.ylabel(r'rate_clusters', fontsize=20)

if savefigures == True:
    plt.savefig(savefigures_directory + 'Power_law_r1_r2_sigma_r/' + model_name + '_Nc_Np_stats.pdf')
else:
    plt.show()
plt.close()
'''





'''
##### To re-make each of the above grid plots as single panels for 2nd year paper:

Model_stats_fields = [r'$f_{\rm sim}/f_{\rm Kepler}$', 'Planet Multiplicities', r'$P$', r'$P_{i+1}/P_i$', r'$\delta$', r'$\/xi$'] #names of the fields associated with the 'Model_stats_...' arrays

#To plot the 2D arrays of model fit statistics on the grid of 'power_law_P' vs. 'power_law_r':
power_law_P_range = np.linspace(-0.5, 0.5, 11)
power_law_r_range = np.linspace(-3., -1.2, 11)

d_y = np.diff(power_law_P_range)[0]
d_x = np.diff(power_law_r_range)[0]
xmin, xmax, ymin, ymax = power_law_r_range[0] - d_x/2., power_law_r_range[-1] + d_x/2., power_law_P_range[0] - d_y/2., power_law_P_range[-1] + d_y/2.
xy_ratio = (xmax - xmin)/(ymax - ymin)

for i,name in enumerate(Model_stats_fields):
    fig = plt.figure(figsize=(5,4))
    plot = GridSpec(1,1,left=0.05,bottom=0.15,right=0.95,top=0.9,wspace=0.1,hspace=0.1)
    ax = plt.subplot(plot[0,0])
    plt.title(name, fontsize=20)
    data_array = Model_stats_P_r_grid[i]
    if i==0:
        data_array = (data_array/Nstars_sim)/(Nplanets_Kepler/Nstars_Kepler)
    plt.imshow(data_array, aspect=xy_ratio, origin='lower', extent=(xmin, xmax, ymin, ymax))
    plt.colorbar()
    plt.xlabel(r'$\alpha_R$', fontsize=20)
    plt.ylabel(r'$\alpha_P$', fontsize=20)
    if savefigures == True:
        plt.savefig(savefigures_directory + 'Power_law_r/Paper_Figures/KS_grids/' + model_name + '_P_r_stats_%s.pdf' % i)
        plt.close()
plt.show()

#To plot the 2D arrays of model fit statistics on the grid of 'num_mutual_hill_radii' vs. 'sigma_logperiod_per_pl_in_cluster':
num_mutual_hill_radii_range = np.linspace(10., 20., 11)
sigma_logperiod_per_pl_in_cluster_range = np.linspace(0.1, 0.4, 11)

d_y = np.diff(num_mutual_hill_radii_range)[0]
d_x = np.diff(sigma_logperiod_per_pl_in_cluster_range)[0]
xmin, xmax, ymin, ymax = sigma_logperiod_per_pl_in_cluster_range[0] - d_x/2., sigma_logperiod_per_pl_in_cluster_range[-1] + d_x/2., num_mutual_hill_radii_range[0] - d_y/2., num_mutual_hill_radii_range[-1] + d_y/2.
xy_ratio = (xmax - xmin)/(ymax - ymin)

for i,name in enumerate(Model_stats_fields):
    fig = plt.figure(figsize=(5,4))
    plot = GridSpec(1,1,left=0.05,bottom=0.15,right=0.95,top=0.9,wspace=0.1,hspace=0.1)
    ax = plt.subplot(plot[0,0])
    plt.title(name, fontsize=20)
    data_array = Model_stats_Nhill_sigma_grid[i]
    if i==0:
        data_array = (data_array/Nstars_sim)/(Nplanets_Kepler/Nstars_Kepler)
    plt.imshow(data_array, aspect=xy_ratio, origin='lower', extent=(xmin, xmax, ymin, ymax))
    plt.colorbar()
    plt.xlabel(r'$\sigma_N$', fontsize=20)
    plt.ylabel(r'$\Delta_c$', fontsize=20)
    if savefigures == True:
        plt.savefig(savefigures_directory + 'Power_law_r/Paper_Figures/KS_grids/' + model_name + '_Nhill_sigma_stats_%s.pdf' % i)
        plt.close()
plt.show()

#To plot the 2D arrays of model fit statistics on the grid of 'sigma_incl' vs. 'sigma_hk':
sigma_incl_range = np.linspace(0., 3., 11)
sigma_hk_range = np.linspace(0., 0.1, 11)

d_y = np.diff(sigma_incl_range)[0]
d_x = np.diff(sigma_hk_range)[0]
xmin, xmax, ymin, ymax = sigma_hk_range[0] - d_x/2., sigma_hk_range[-1] + d_x/2., sigma_incl_range[0] - d_y/2., sigma_incl_range[-1] + d_y/2.
xy_ratio = (xmax - xmin)/(ymax - ymin)

for i,name in enumerate(Model_stats_fields):
    fig = plt.figure(figsize=(5,4))
    plot = GridSpec(1,1,left=0.05,bottom=0.15,right=0.95,top=0.9,wspace=0.1,hspace=0.1)
    ax = plt.subplot(plot[0,0])
    plt.title(name, fontsize=20)
    data_array = Model_stats_sigma_incl_sigma_hk_grid[i]
    if i==0:
        data_array = (data_array/Nstars_sim)/(Nplanets_Kepler/Nstars_Kepler)
    plt.imshow(data_array, aspect=xy_ratio, origin='lower', extent=(xmin, xmax, ymin, ymax))
    plt.colorbar()
    plt.xlabel(r'$\sigma_e$', fontsize=20)
    plt.ylabel(r'$\sigma_i$', fontsize=20)
    if savefigures == True:
        plt.savefig(savefigures_directory + 'Power_law_r/Paper_Figures/KS_grids/' + model_name + '_sigma_incl_sigma_hk_stats_%s.pdf' % i)
        plt.close()
plt.show()

#To plot the 2D arrays of model fit statistics on the grid of 'max_radius' vs. 'mr_power_law':
max_radius_range = np.linspace(2., 12., 11)
mr_range = np.linspace(1.5, 3.5, 11)

d_y = np.diff(max_radius_range)[0]
d_x = np.diff(mr_range)[0]
xmin, xmax, ymin, ymax = mr_range[0] - d_x/2., mr_range[-1] + d_x/2., max_radius_range[0] - d_y/2., max_radius_range[-1] + d_y/2.
xy_ratio = (xmax - xmin)/(ymax - ymin)

for i,name in enumerate(Model_stats_fields):
    fig = plt.figure(figsize=(5,4))
    plot = GridSpec(1,1,left=0.05,bottom=0.15,right=0.95,top=0.9,wspace=0.1,hspace=0.1)
    ax = plt.subplot(plot[0,0])
    plt.title(name, fontsize=20)
    data_array = Model_stats_max_radius_mr_grid[i]
    if i==0:
        data_array = (data_array/Nstars_sim)/(Nplanets_Kepler/Nstars_Kepler)
    plt.imshow(data_array, aspect=xy_ratio, origin='lower', extent=(xmin, xmax, ymin, ymax))
    plt.colorbar()
    plt.xlabel(r'$\alpha_{mr}$', fontsize=20)
    plt.ylabel(r'$R_{p,\rm max}$ $(R_\odot)$', fontsize=20)
    if savefigures == True:
        plt.savefig(savefigures_directory + 'Power_law_r/Paper_Figures/KS_grids/' + model_name + '_max_radius_mr_stats_%s.pdf' % i)
        plt.close()
plt.show()

#To plot the 2D arrays of model fit statistics on the grid of 'log_rate_clusters' vs. 'log_rate_planets_per_cluster':
Model_stats_fields = [r'$f_{\rm sim}/f_{\rm Kepler}$', r'$N_{\rm obs,singles}$', 'Planet Multiplicities', r'$P$', r'$P_{i+1}/P_i$', r'$\delta$', r'$\/xi$'] #names of the fields associated with the 'Model_stats_...' arrays

Nc_range = np.linspace(1., 3., 11)
Np_range = np.linspace(1., 3., 11)

d_y = np.diff(Nc_range)[0]
d_x = np.diff(Np_range)[0]
xmin, xmax, ymin, ymax = Np_range[0] - d_x/2., Np_range[-1] + d_x/2., Nc_range[0] - d_y/2., Nc_range[-1] + d_y/2.
xy_ratio = (xmax - xmin)/(ymax - ymin)

for i,name in enumerate(Model_stats_fields):
    fig = plt.figure(figsize=(5,4))
    plot = GridSpec(1,1,left=0.05,bottom=0.15,right=0.95,top=0.9,wspace=0.1,hspace=0.1)
    ax = plt.subplot(plot[0,0])
    plt.title(name, fontsize=20)
    data_array = Model_stats_Nc_Np_grid[i]
    if i==0:
        data_array = (data_array/Nstars_sim)/(Nplanets_Kepler/Nstars_Kepler)
    plt.imshow(data_array, aspect=xy_ratio, origin='lower', extent=(xmin, xmax, ymin, ymax))
    plt.colorbar()
    plt.xlabel(r'$\lambda_p$', fontsize=20)
    plt.ylabel(r'$\lambda_c$', fontsize=20)
    if savefigures == True:
        plt.savefig(savefigures_directory + 'Power_law_r/Paper_Figures/KS_grids/' + model_name + '_Nc_Np_stats_%s.pdf' % i)
        plt.close()
plt.show()
'''





'''
#To plot the 2D arrays of model fit statistics on the grid of 'power_law_r1' vs. 'power_law_r2' for broken power law model of radii:
Model_stats_fields = ['Number of observed planets', 'K-S multiplicities', 'K-S periods', 'K-S period ratios', 'K-S transit durations', 'K-S xi', 'K-S transit depths', 'K-S transit depth ratios'] #names of the fields associated with the 'Model_stats_...' arrays

r1_range = np.linspace(-6., 3., 11)
r2_range = np.linspace(-6., 3., 11)

d_y = np.diff(r1_range)[0]
d_x = np.diff(r2_range)[0]
xmin, xmax, ymin, ymax = r2_range[0] - d_x/2., r2_range[-1] + d_x/2., r1_range[0] - d_y/2., r1_range[-1] + d_y/2.
xy_ratio = (xmax - xmin)/(ymax - ymin)

fig = plt.figure(figsize=(16,8))
plot = GridSpec(2,4,left=0.075,bottom=0.115,right=0.95,top=0.925,wspace=0.25,hspace=0.4)

plot_rows = [0,0,0,0,1,1,1,1]
plot_cols = [0,1,2,3,0,1,2,3]
for i,name in enumerate(Model_stats_fields):
    ax = plt.subplot(plot[plot_rows[i],plot_cols[i]])
    plt.title(name, fontsize=20)
    plt.imshow(Model_stats_r1_r2_grid[i], aspect=xy_ratio, origin='lower', extent=(xmin, xmax, ymin, ymax))
    plt.colorbar()
    plt.xlabel(r'power_law_r2', fontsize=20)
    plt.ylabel(r'power_law_r1', fontsize=20)

if savefigures == True:
    plt.savefig(savefigures_directory + 'Power_law_r1_r2_sigma_r/' + model_name + '_r1_r2_rbreak2_stats.pdf')
else:
    plt.show()
plt.close()
'''
    
'''
#To plot the 2D arrays of model fit statistics on the grid of 'power_law_r1' vs. 'power_law_r2' for broken power law model of radii:
r1_range = np.linspace(-3., 3., 11)
r2_range = np.linspace(-3., 3., 11)

d_y = np.diff(r1_range)[0]
d_x = np.diff(r2_range)[0]
xmin, xmax, ymin, ymax = r2_range[0] - d_x/2., r2_range[-1] + d_x/2., r1_range[0] - d_y/2., r1_range[-1] + d_y/2.
xy_ratio = (xmax - xmin)/(ymax - ymin)

fig = plt.figure(figsize=(16,8))
plot = GridSpec(2,4,left=0.075,bottom=0.115,right=0.95,top=0.925,wspace=0.25,hspace=0.4)

plot_rows = [0,0,0,0,1,1,1]
plot_cols = [0,1,2,3,0,1,2]
for i,name in enumerate(Model_stats_fields):
    ax = plt.subplot(plot[plot_rows[i],plot_cols[i]])
    plt.title(name, fontsize=20)
    plt.imshow(Model_stats_r1_r2_rbreak3_grid[i], aspect=xy_ratio, origin='lower', extent=(xmin, xmax, ymin, ymax))
    plt.colorbar()
    plt.xlabel(r'power_law_r2', fontsize=20)
    plt.ylabel(r'power_law_r1', fontsize=20)

if savefigures == True:
    plt.savefig(savefigures_directory + 'Power_law_r1_r2/' + model_name + '_r1_r2_rbreak3_stats.pdf')
else:
    plt.show()
plt.close()
'''





'''
#To plot the 2D arrays of model fit statistics on the grid of 'power_law_r' vs. 'sigma_log_radius_in_cluster' for clustered (lognormal) radii, but still a single power law:
Model_stats_fields = ['Number of observed planets', 'Number of observed singles', 'K-S multiplicities', 'K-S periods', 'K-S period ratios', 'K-S transit depths', 'K-S xi'] #names of the fields associated with the 'Model_stats_...' arrays

r_range = np.linspace(-3., -1.2, 11)
sigma_r_range = np.linspace(0.05, 0.5, 11)

d_y = np.diff(r_range)[0]
d_x = np.diff(sigma_r_range)[0]
xmin, xmax, ymin, ymax = sigma_r_range[0] - d_x/2., sigma_r_range[-1] + d_x/2., r_range[0] - d_y/2., r_range[-1] + d_y/2.
xy_ratio = (xmax - xmin)/(ymax - ymin)

fig = plt.figure(figsize=(16,8))
plot = GridSpec(2,4,left=0.075,bottom=0.115,right=0.95,top=0.925,wspace=0.25,hspace=0.4)

plot_rows = [0,0,0,0,1,1,1]
plot_cols = [0,1,2,3,0,1,2]
for i,name in enumerate(Model_stats_fields):
    ax = plt.subplot(plot[plot_rows[i],plot_cols[i]])
    plt.title(name, fontsize=20)
    plt.imshow(Model_stats_r_sigma_r_grid[i], aspect=xy_ratio, origin='lower', extent=(xmin, xmax, ymin, ymax))
    plt.colorbar()
    plt.xlabel(r'sigma_log_radius_in_cluster', fontsize=20)
    plt.ylabel(r'power_law_r', fontsize=20)

if savefigures == True:
    plt.savefig(savefigures_directory + 'Power_law_r_sigma_r/' + model_name + '_r_sigma_r_stats.pdf')
else:
    plt.show()
plt.close()
'''





#'''
##### To re-make each of the above grid plots as single panels for Comprehensive exam paper:

Model_stats_fields = [r'$f_{\rm sim}/f_{\rm Kepler}$', 'Planet Multiplicities', r'$P$', r'$P_{i+1}/P_i$', r'$t_{\rm dur}$', r'$\xi$', r'$\delta$', r'$\delta_{i+1}/\delta_i$'] #names of the fields associated with the 'Model_stats_...' arrays

#To plot the 2D arrays of model fit statistics on the grid of 'power_law_r1' vs. 'power_law_r2':
r1_range = np.linspace(-6., 3., 11)
r2_range = np.linspace(-6., 3., 11)

d_y = np.diff(r1_range)[0]
d_x = np.diff(r2_range)[0]
xmin, xmax, ymin, ymax = r2_range[0] - d_x/2., r2_range[-1] + d_x/2., r1_range[0] - d_y/2., r1_range[-1] + d_y/2.
xy_ratio = (xmax - xmin)/(ymax - ymin)

for i,name in enumerate(Model_stats_fields):
    fig = plt.figure(figsize=(5,4))
    plot = GridSpec(1,1,left=0.05,bottom=0.15,right=0.95,top=0.9,wspace=0.1,hspace=0.1)
    ax = plt.subplot(plot[0,0])
    plt.title(name, fontsize=20)
    data_array = Model_stats_r1_r2_grid[i]
    if i==0:
        data_array = (data_array/(Nstars_sim/cos_factor))/(Nplanets_Kepler/Nstars_Kepler)
        plt.imshow(data_array, aspect=xy_ratio, vmin=0.5, vmax=2., origin='lower', extent=(xmin, xmax, ymin, ymax))
    else:
        plt.imshow(data_array, aspect=xy_ratio, vmin=0., vmax=0.2, origin='lower', extent=(xmin, xmax, ymin, ymax))
    plt.colorbar()
    plt.xlabel(r'$\alpha_{R2}$', fontsize=20)
    plt.ylabel(r'$\alpha_{R1}$', fontsize=20)
    if savefigures == True:
        plt.savefig(savefigures_directory + 'Power_law_r1_r2_sigma_r/Paper_Figures/KS_grids/' + model_name + '_r1_r2_rbreak2_stats_%s.pdf' % i)
        plt.close()
plt.show()

#To plot the 2D arrays of model fit statistics on the grid of 'num_mutual_hill_radii' vs. 'sigma_logperiod_per_pl_in_cluster':
num_mutual_hill_radii_range = np.linspace(5., 20., 11)
sigma_logperiod_per_pl_in_cluster_range = np.linspace(0.05, 0.5, 11)

d_y = np.diff(num_mutual_hill_radii_range)[0]
d_x = np.diff(sigma_logperiod_per_pl_in_cluster_range)[0]
xmin, xmax, ymin, ymax = sigma_logperiod_per_pl_in_cluster_range[0] - d_x/2., sigma_logperiod_per_pl_in_cluster_range[-1] + d_x/2., num_mutual_hill_radii_range[0] - d_y/2., num_mutual_hill_radii_range[-1] + d_y/2.
xy_ratio = (xmax - xmin)/(ymax - ymin)

for i,name in enumerate(Model_stats_fields):
    fig = plt.figure(figsize=(5,4))
    plot = GridSpec(1,1,left=0.05,bottom=0.15,right=0.95,top=0.9,wspace=0.1,hspace=0.1)
    ax = plt.subplot(plot[0,0])
    plt.title(name, fontsize=20)
    data_array = Model_stats_Nhill_sigma_grid[i]
    if i==0:
        data_array = (data_array/(Nstars_sim/cos_factor))/(Nplanets_Kepler/Nstars_Kepler)
        plt.imshow(data_array, aspect=xy_ratio, vmin=0.5, vmax=2., origin='lower', extent=(xmin, xmax, ymin, ymax))
    else:
        plt.imshow(data_array, aspect=xy_ratio, vmin=0., vmax=0.2, origin='lower', extent=(xmin, xmax, ymin, ymax))
    plt.colorbar()
    plt.xlabel(r'$\sigma_N$', fontsize=20)
    plt.ylabel(r'$\Delta_c$', fontsize=20)
    if savefigures == True:
        plt.savefig(savefigures_directory + 'Power_law_r1_r2_sigma_r/Paper_Figures/KS_grids/' + model_name + '_Nhill_sigma_stats_%s.pdf' % i)
        plt.close()
plt.show()

#To plot the 2D arrays of model fit statistics on the grid of 'sigma_incl' vs. 'sigma_hk':
sigma_incl_range = np.linspace(0., 5., 11)
sigma_hk_range = np.linspace(0., 0.2, 11)

d_y = np.diff(sigma_incl_range)[0]
d_x = np.diff(sigma_hk_range)[0]
xmin, xmax, ymin, ymax = sigma_hk_range[0] - d_x/2., sigma_hk_range[-1] + d_x/2., sigma_incl_range[0] - d_y/2., sigma_incl_range[-1] + d_y/2.
xy_ratio = (xmax - xmin)/(ymax - ymin)

for i,name in enumerate(Model_stats_fields):
    fig = plt.figure(figsize=(5,4))
    plot = GridSpec(1,1,left=0.05,bottom=0.15,right=0.95,top=0.9,wspace=0.1,hspace=0.1)
    ax = plt.subplot(plot[0,0])
    plt.title(name, fontsize=20)
    data_array = Model_stats_sigma_incl_sigma_hk_grid[i]
    if i==0:
        data_array = (data_array/(Nstars_sim/cos_factor))/(Nplanets_Kepler/Nstars_Kepler)
        plt.imshow(data_array, aspect=xy_ratio, vmin=0.5, vmax=2., origin='lower', extent=(xmin, xmax, ymin, ymax))
    else:
        plt.imshow(data_array, aspect=xy_ratio, vmin=0., vmax=0.2, origin='lower', extent=(xmin, xmax, ymin, ymax))
    plt.colorbar()
    plt.xlabel(r'$\sigma_e$', fontsize=20)
    plt.ylabel(r'$\sigma_i$', fontsize=20)
    if savefigures == True:
        plt.savefig(savefigures_directory + 'Power_law_r1_r2_sigma_r/Paper_Figures/KS_grids/' + model_name + '_sigma_incl_sigma_hk_stats_%s.pdf' % i)
        plt.close()
plt.show()

#To plot the 2D arrays of model fit statistics on the grid of 'mr_power_law' vs. 'sigma_log_radius_in_cluster':
mr_range = np.linspace(1.5, 3.5, 11)
sigmaR_range = np.linspace(0.1, 1.0, 11)

d_y = np.diff(mr_range)[0]
d_x = np.diff(sigmaR_range)[0]
xmin, xmax, ymin, ymax = sigmaR_range[0] - d_x/2., sigmaR_range[-1] + d_x/2., mr_range[0] - d_y/2., mr_range[-1] + d_y/2.
xy_ratio = (xmax - xmin)/(ymax - ymin)

for i,name in enumerate(Model_stats_fields):
    fig = plt.figure(figsize=(5,4))
    plot = GridSpec(1,1,left=0.05,bottom=0.15,right=0.95,top=0.9,wspace=0.1,hspace=0.1)
    ax = plt.subplot(plot[0,0])
    plt.title(name, fontsize=20)
    data_array = Model_stats_mr_sigmaR_grid[i]
    if i==0:
        data_array = (data_array/(Nstars_sim/cos_factor))/(Nplanets_Kepler/Nstars_Kepler)
        plt.imshow(data_array, aspect=xy_ratio, vmin=0.5, vmax=2., origin='lower', extent=(xmin, xmax, ymin, ymax))
    else:
        plt.imshow(data_array, aspect=xy_ratio, vmin=0., vmax=0.2, origin='lower', extent=(xmin, xmax, ymin, ymax))
    plt.colorbar()
    plt.xlabel(r'$\sigma_R$', fontsize=20)
    plt.ylabel(r'$\alpha_{mr}$', fontsize=20)
    if savefigures == True:
        plt.savefig(savefigures_directory + 'Power_law_r1_r2_sigma_r/Paper_Figures/KS_grids/' + model_name + '_mr_sigmaR_stats_%s.pdf' % i)
        plt.close()
plt.show()

#To plot the 2D arrays of model fit statistics on the grid of 'log_rate_clusters' vs. 'log_rate_planets_per_cluster':
Model_stats_fields = [r'$f_{\rm sim}/f_{\rm Kepler}$', r'$N_{p,\rm obs,singles}$', 'Planet Multiplicities', r'$P$', r'$P_{i+1}/P_i$', r'$t_{\rm dur}$', r'$\xi$', r'$\delta$', r'$\delta_{i+1}/\delta_i$'] #names of the fields associated with the 'Model_stats_...' arrays

Nc_range = np.linspace(1., 5., 11)
Np_range = np.linspace(1., 5., 11)

d_y = np.diff(Nc_range)[0]
d_x = np.diff(Np_range)[0]
xmin, xmax, ymin, ymax = Np_range[0] - d_x/2., Np_range[-1] + d_x/2., Nc_range[0] - d_y/2., Nc_range[-1] + d_y/2.
xy_ratio = (xmax - xmin)/(ymax - ymin)

for i,name in enumerate(Model_stats_fields):
    fig = plt.figure(figsize=(5,4))
    plot = GridSpec(1,1,left=0.05,bottom=0.15,right=0.95,top=0.9,wspace=0.1,hspace=0.1)
    ax = plt.subplot(plot[0,0])
    plt.title(name, fontsize=20)
    data_array = Model_stats_Nc_Np_grid[i]
    if i==0:
        data_array = (data_array/(Nstars_sim/cos_factor))/(Nplanets_Kepler/Nstars_Kepler)
        plt.imshow(data_array, aspect=xy_ratio, vmin=0.5, vmax=2., origin='lower', extent=(xmin, xmax, ymin, ymax))
    else:
        plt.imshow(data_array, aspect=xy_ratio, vmin=0., vmax=0.2, origin='lower', extent=(xmin, xmax, ymin, ymax))
    plt.colorbar()
    plt.xlabel(r'$\lambda_p$', fontsize=20)
    plt.ylabel(r'$\lambda_c$', fontsize=20)
    if savefigures == True:
        plt.savefig(savefigures_directory + 'Power_law_r1_r2_sigma_r/Paper_Figures/KS_grids/' + model_name + '_Nc_Np_stats_%s.pdf' % i)
        plt.close()
plt.show()
#'''
