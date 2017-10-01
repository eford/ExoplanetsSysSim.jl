## ExoplanetsSysSim/examples/hsu_et_al/sim_param_plots.jl
## (c) 2017 Keir A. Ashby
# Takes any number of SimParams and, calculates observed and physical catalogs with summary statistics and returns:
#   1.A period-radius distribution of the observed, physical, and true catalogs for each sim_param individually 
#   2.A period-radius distribution of the observed and true catalogs for each sim_param individually 
#   3.An overlaid comparison period-radius distribution of the observed and true catalogs for every possible combination of two input sim_params
#   4.An overlaid comparison period-radius distribution of the observed and true catalogs for every input sim_param
#   5.The percentage of observed planets to true planets in each period-radius bin
#   6.CDF plots for radii and periods

using PyPlot
using ExoplanetsSysSim
using Combinatorics

include(joinpath(Pkg.dir(),"ExoplanetsSysSim","src","constants.jl"))
include(joinpath(Pkg.dir(),"ExoplanetsSysSim","examples","hsu_etal_2017", "abc_setup_christiansen.jl"))
include(joinpath(Pkg.dir(),"ExoplanetsSysSim","examples","hsu_etal_2017", "christiansen_func.jl"))
include(joinpath(Pkg.dir(),"ExoplanetsSysSim","examples","hsu_etal_2017", "calc_sum_stats_pr.jl"))

#see general file definition at the top for a description of this function's purpose.
function sim_param_plots(have_sim_params::Bool; sim_arr::Array{Any,1} = [], num_param::Int64 = 0, print_rates=false) 	#if the user already has SimParams to put in, the first input will be true and the second will be an array (even if there's only one) of your SimParams. If you don't already have one, put false as your first input, and how many SimParams you'd like to generate. If you want to see the ratios of observed to kepler planets in all the bins, set print_rates=true) 
  if have_sim_params == false 						#we will make SimParams for the user if she doesn't already have one
    sim_arr = setup_sim_arr(num_param) 					#returns an array of SimParams
  end
  sim_true = setup_sim_arr(1)[1] 					#generates a SimParam that will be used to get the true Kepler catalog.
  ss_arr, ss_true = setup_sum_stats(sim_arr, sim_true) 			#generates physical and observed catalogs and returns summary statistics for each SimParam and for the true Kepler catalog 
  num_cats = length(ss_arr)
  n=1									#iterative integer to keep track of figures
  pr_lists = make_pr_lists(num_cats, ss_arr, ss_true)			#Makes lists of all the periods and radii for each catalog
  n = all_cats_pr(n, num_cats, pr_lists)				#A period-radius distribution of the observed, physical, and true catalogs for each sim_param individually 
  n = all_cats_pr(n, num_cats, pr_lists, false)				#A period-radius distribution of the observed and true catalogs for each sim_param individually 
  if num_cats > 1
    n = compare_2_obs_pr(n, num_cats, pr_lists)				#An overlaid comparison period-radius distribution of the observed and true catalogs for every possible combination of two input sim_params
  end
  if num_cats > 2
    n = compare_all_obs_pr(n, num_cats, pr_lists)			#An overlaid comparison period-radius distribution of the observed and true catalogs for every input sim_param
  end
  ratios = pr_rates(pr_lists["periods_obs_1"], pr_lists["radii_obs_1"], pr_lists["periods_true"], pr_lists["radii_true"], pr_lists["p_lim_arr"], pr_lists["r_lim_arr"]) #The percentage of observed planets to true planets in each period-radius bin  
  n = cdf_plot(n, pr_lists["periods_obs_1"], pr_lists["periods_true"], pr_lists["p_lim_arr"], "Periods") 	#CDF plots for periods
  n = cdf_plot(n, pr_lists["radii_obs_1"], pr_lists["radii_true"], pr_lists["r_lim_arr"], "Radii")		#CDF plots for radii
  return(ratios)							#Normally, pr_rates would just be called last so it printed ratios at the end automatically. However, cdf_plot for radii alters the radii_obs_1 list so that it doesn't work in pr_rates. I'm working on how to alter the list in cdf_plot without altering it globally.
end

function setup_sim_arr(num_param::Int64 = 1) 				#takes an integer and returns an array of that many Sim Params 
  sim_arr = Any[]
  for i in 1:num_param
    sim_param = setup_sim_param_christiansen() 				#sim-param used in christiansen_func.jl
    hsu_final_rates = [0.039  0.4    2.6   5.1   5.4  3.9  1.6  56  87 
		       0.1    0.3    1     3.2   3.6  7.9  1.7  5.5 21
		       0.088  0.18   0.52  1.5   2.8  3.1  4.1  3.2 8.5
		       0.053  0.13   0.63  1.6   2.4  2.4  2.3  4.4 7.6
		       0.047  0.13   0.74  1.1   1.9  1.9  2.3  2.7 2.2
		       0.036  0.057  0.38  0.8   1.6  1.2  2    2.4 3.5
		       0.0096 0.059  0.29  1.1   2.5  3.7  3.5  3.6 6.4
   		       0.007  0.042  0.22  0.93  2.3  2.8  3.5  3.1 3.2
		       0.0044 0.029  0.14  0.6   1.4  1.9  2.3  3.4 5
		       0.0086 0.012  0.13  0.23  0.48 0.82 0.91 1.5 1.8
		       0.0033 0.0068 0.057 0.1   0.19 0.23 0.76 0.7 1.9
       		       0.0045 0.012  0.082 0.16  0.23 0.41 0.32 1.2 2.4
		       0.0045 0.031  0.12  0.075 0.16 0.18 0.33 0.5 0.59]*0.01	#occurence rates over all bins drawn from results in Hsu et al 2017
    hsu_radii = [0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 4, 6, 8, 12, 16]	#limits for the radius bins
    hsu_periods = [0.5, 1.25, 2.5, 5, 10, 20, 40, 80, 160, 320]			#limits for the period bins
    update_param(sim_param, "obs_par", hsu_final_rates)
    update_param(sim_param, "r_lim_arr", hsu_radii*ExoplanetsSysSim.earth_radius)
    update_param(sim_param, "p_lim_arr", hsu_periods)
    update_param(sim_param, "num_targets_sim_pass_one", 150000)
    update_param(sim_param, "num_kepler_targets", 150000)
    push!(sim_arr,sim_param)
  end
  return(sim_arr)
end

#takes an array of Sim Params and returns summary statistics for each. It will also take one Sim Param and make summary statistics for the true Kepler Catalog
function setup_sum_stats(sim_arr::Array{Any,1}, sim_true::ExoplanetsSysSim.SimParam)
  setup_MES()
  cat_true = setup_actual_planet_candidate_catalog(setup_star_table_christiansen(sim_true), sim_true)
  ss_true = calc_sum_stats_pr(cat_true, sim_true, true)
  ss_arr = Array{Any}(length(sim_arr))
  for i in 1:length(sim_arr)
    cat_phys = generate_kepler_physical_catalog(sim_arr[i])
    cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys, sim_arr[i])  
    cat_obs = ExoplanetsSysSim.observe_kepler_targets_single_obs(cat_phys_cut, sim_arr[i])
    ss_arr[i] = calc_sum_stats_pr(cat_obs, sim_arr[i], false, cat_phys)
  end
  return(ss_arr, ss_true)
end

function make_pr_lists(num_cats, ss_arr::Array{Any,1}, ss_true::ExoplanetsSysSim.CatalogSummaryStatistics)	#Makes lists of all the periods and radii for each catalog
  p_lim_arr = [0.5, 1.25, 2.5, 5, 10, 20, 40, 80, 160, 320]				#bin boundaries for periods
  r_lim_arr = [0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 4, 6, 8, 12, 16, 20]		#bin boundaries for radii
  pr_lists = Dict("p_lim_arr"=> p_lim_arr, "r_lim_arr"=> r_lim_arr)
  for i in 1:num_cats
    pr_lists["periods_obs_$i"] = ss_arr[i].stat["period_obs_list"] 			#all the periods in this observable catalog
    pr_lists["radii_obs_$i"] = ss_arr[i].stat["radius_obs_list"]./earth_radius 		#all the radii in this observable catalog
    pr_lists["periods_phys_$i"] = ss_arr[i].stat["period_phys_list"] 			#all the periods in this physical catalog
    pr_lists["radii_phys_$i"] = ss_arr[i].stat["radius_phys_list"]./earth_radius 	#all the radii in this physical catalog 
    if length(pr_lists["periods_phys_$i"])==0 || length(pr_lists["periods_obs_$i"])==0 	#test to make sure we have something to plot
      println("One of your catalogs is empty")
      return()
    end  
  end
  pr_lists["periods_true"] = ss_true.stat["period_list"] 				#all the periods in this true catalog
  pr_lists["radii_true"] = ss_true.stat["radius_list"]./earth_radius 			#all the radii in this true catalog
  return(pr_lists)
end

#A period-radius distribution of the observed, physical, and true catalogs for each sim_param individually 
function all_cats_pr(n::Int64, num_cats::Int64, pr_lists::Dict{String,Array{Float64,1}}, plot_phys = true)
  for i in 1:num_cats
    figure(n)
    if plot_phys == true
      loglog(pr_lists["periods_phys_$i"], pr_lists["radii_phys_$i"], ".", ms = 0.5, color = "orange", label = "Physical")
    end
    loglog(pr_lists["periods_obs_$i"], pr_lists["radii_obs_$i"], ".", ms = 3, label = "Simulated Observed")
    loglog(pr_lists["periods_true"], pr_lists["radii_true"], "x", ms = 2, color = "purple", label = "Kepler Observed")
    title("Periods vs. Radii $i", fontsize = 18)
    axis([0.0, 350.0, 0.4, 22.0])								
    xticks(pr_lists["p_lim_arr"], pr_lists["p_lim_arr"])
    yticks(pr_lists["r_lim_arr"], pr_lists["r_lim_arr"])
    xlabel("Periods (in days)", fontsize=14)
    ylabel("Radius (in Earth radii)", fontsize = 14)
    legend(loc = "lower right", fontsize = 10)
    if plot_phys == true 
     savefig("all_cats_pr_$i.pdf")
    else
     savefig("obs_cats_pr_$i.pdf")
    end
    n+=1
  end
  return(n)
end

#An overlaid comparison period-radius distribution of the observed and true catalogs for every possible combination of two input sim_params
function compare_2_obs_pr(n::Int64, num_cats::Int64, pr_lists::Dict{String,Array{Float64,1}})
  for (a,b) in combinations(collect(1:num_cats),2)		
    figure(n)
    loglog(pr_lists["periods_true"], pr_lists["radii_true"], "x", ms = 2, color = "purple", label = "Kepler Observed")
    loglog(pr_lists["periods_obs_$a"], pr_lists["radii_obs_$a"], ".", ms = 3, label = "Simulated Observed: $a")
    loglog(pr_lists["periods_obs_$b"], pr_lists["radii_obs_$b"], ".", ms = 3, label = "Simulated Observed: $b")
    title("Periods vs. Radii Comparison: $a & $b", fontsize = 18)
    axis([0.0, 350.0, 0.4, 22.0])								
    xticks(pr_lists["p_lim_arr"], pr_lists["p_lim_arr"])
    yticks(pr_lists["r_lim_arr"], pr_lists["r_lim_arr"])
    xlabel("Periods (in days)", fontsize = 14)
    ylabel("Radius (in Earth radii)", fontsize = 14)
    legend(loc = "lower right", fontsize = 10)
    n+=1
#    savefig("compare_2_obs_pr_$a_&_$b.pdf")
  end
  return(n)
end

#An overlaid comparison period-radius distribution of the observed and true catalogs for every possible combination of two input sim_params
function compare_all_obs_pr(n::Int64, num_cats::Int64, pr_lists::Dict{String,Array{Float64,1}}) 	
  figure(n)
  loglog(pr_lists["periods_true"], pr_lists["radii_true"], "x", ms = 2, color = "purple", label = "Kepler Observed")
  for i in 1:num_cats					
    loglog(pr_lists["periods_obs_$i"], pr_lists["radii_obs_$i"], ".", ms = 3, label = "Simulated Observed: $i")
  end
  title("Periods vs. Radii Comparison: All", fontsize = 18)
  axis([0.0, 350.0, 0.4, 22.0])								
  xticks(pr_lists["p_lim_arr"], pr_lists["p_lim_arr"])
  yticks(pr_lists["r_lim_arr"], pr_lists["r_lim_arr"])
  xlabel("Periods (in days)", fontsize=14)
  ylabel("Radius (in Earth radii)", fontsize = 14)
  legend(loc = "lower right", fontsize = 10)
#  savefig("compare_all_obs_pr.pdf")
  n+=1
  return(n)
end

#CDF plots for radii and periods
function cdf_plot(n::Int64, obs_list::Array{Float64,1}, true_list::Array{Float64,1}, lim_arr::Array{Float64,1}, cdf_type::String) 
  figure(n)
  obs_list = trim_catalog(obs_list, lim_arr[1], lim_arr[length(lim_arr)])	#we want to restrict our catalog to stay within our bin ranges
  true_list = trim_catalog(true_list, lim_arr[1], lim_arr[length(lim_arr)])	#we want to restrict our catalog to stay within our bin ranges
  semilogx(obs_list, collect(linspace(0.0,1.0,length(obs_list))), label = "Simulated Observed")  
  semilogx(true_list, collect(linspace(0.0,1.0,length(true_list))), label = "Kepler Observed")  
  title("Cumulative Distribution of $cdf_type")
  xlabel("Periods (in days)")
  ylabel("Probability")
  legend(loc = "lower right", fontsize = 10)
#  savefig("cdf_$cdf_type.pdf")
  n+=1
  return(n)
end

#This is a function for testing the accuracy of our model for how many planets Kepler would see. Returns a 2X2 array  of observed/true rates in fraction or decimal form
function pr_rates(periods_obs::Array{Float64,1}, radii_obs::Array{Float64,1}, periods_true::Array{Float64,1}, radii_true::Array{Float64,1}, p_lim_arr::Array{Float64,1}, r_lim_arr::Array{Float64,1}, percent::Bool = true) #if percent is false, the rates will be printed as fractions
  rates_obs = zeros(length(p_lim_arr)-1,length(r_lim_arr)-1)			#table to hold the number of observed planets in each bin
  rates_true = zeros(length(p_lim_arr)-1,length(r_lim_arr)-1)			#table to hold the number of true planets in each bin
  ratios = Array{String}(length(p_lim_arr)-1,length(r_lim_arr)-1)		#table to hold the ratios of observed/true planets in each bin
  for i in 1:length(p_lim_arr)-1						#for each column of periods
    f(p) = p>p_lim_arr[i] && p<p_lim_arr[i+1]    				
    p_bins_obs = find(f,periods_obs)						#this will hold all the indices of planets that fall into the ith period bin
    r_bins_obs = zeros(length(p_bins_obs))					#this will hold the radius corresponding to each exoplanet whose index is held in p_bins_obs
    for j in 1:length(p_bins_obs)
      r_bins_obs[j]=radii_obs[p_bins_obs][j]
    end
    for k in 1:length(r_lim_arr)-1						#for each row of radii in the ith column of periods
      g(r) = r>r_lim_arr[k] && r<r_lim_arr[k+1]
      rates_obs[i,length(r_lim_arr)-k] = convert(Int64,length(find(g,r_bins_obs))) 
    end
  end
#  for h in 1:length(r_lim_arr)-1 						#shows rates_obs 
#    println(rates[:,h])
#  end

  for i in 1:length(p_lim_arr)-1						#for each column of periods
    f(p) = p>p_lim_arr[i] && p<p_lim_arr[i+1]    
    p_bins_true = find(f,periods_true)						#this will hold all the indices of planets that fall into the ith period bin
    r_bins_true = zeros(length(p_bins_true))					#this will hold the radius corresponding to each exoplanet whose index is held in p_bins_true
    for j in 1:length(p_bins_true)
      r_bins_true[j]=radii_true[p_bins_true][j]
    end
    for k in 1:length(r_lim_arr)-1						#for each row of radii in the ith column of periods
      g(r) = r>r_lim_arr[k] && r<r_lim_arr[k+1]
      rates_true[i,length(r_lim_arr)-k] = convert(Int64,length(find(g,r_bins_true)))
    end
  end       
#  for h in 1:length(r_lim_arr)-1						#shows rates_true
#    println(rates_true[:,h])
#  end 	

  for i in 1:length(p_lim_arr)-1
    for j in 1:length(r_lim_arr)-1
      if percent == true
	if rates_obs[i,j] == 0.0 && rates_true[i,j] == 0.0
	  ratios[i,j] = "NaN"
        elseif rates_true[i,j] == 0.0
	  ratios[i,j] = "Inf"
	else
	  ratios[i,j] = string(floor(Int64,rates_obs[i,j]/rates_true[i,j]*100), "%")
	end
      else    									#returns the ratios in fraction form
        ratios[i,j] = string(convert(Int64,rates_obs[i,j]), "/", convert(Int64,rates_true[i,j]))
      end
    end
  end
  return ratios 	
end

function trim_catalog(list::Array{Float64,1},min::Float64,max::Float64)		#removes all values in a list to stay within bin ranges
  sort!(list)
  while list[1] < min
    shift!(list)
  end
  while list[length(list)] > max
    pop!(list)
  end
  return(list)
end

function test_period_radius_plots()
  sim_param_plots(false, num_param = 3)
end
