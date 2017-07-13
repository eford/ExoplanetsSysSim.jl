using ExoplanetsSysSim
using StatsFuns
using JLD
using DataFrames
#import ExoplanetsSysSim.StellarTable.df
#import ExoplanetsSysSim.StellarTable.usable
import Compat: UTF8String, ASCIIString

## simulation_parameters
function setup_sim_param_model(args::Vector{ASCIIString} = Array{ASCIIString}(0) )   # allow this to take a list of parameter (e.g., from command line)
  sim_param = SimParam()
  # How many tatrges to generate
  #add_param_fixed(sim_param,"num_targets_sim_pass_one",150969)                      # Note this is used for the number of stars in the simulations, not necessarily related to number of Kepler targets
  add_param_fixed(sim_param,"num_targets_sim_pass_one",15096)                      # Note this is used for the number of stars in the simulations, not necessarily related to number of Kepler targets
  add_param_fixed(sim_param,"num_kepler_targets",150969)                            # Note this is used for the number of Kepler targets for the observational catalog

  # For generating target star properties
  add_param_fixed(sim_param,"generate_kepler_target",generate_kepler_target_from_table)
  add_param_fixed(sim_param,"star_table_setup",setup_star_table_christiansen)
  add_param_fixed(sim_param,"stellar_catalog","q1_q17_dr25_stellar.jld")
  # add_param_fixed(sim_param,"generate_kepler_target",ExoplanetsSysSim.generate_kepler_target_simple)  # An alternative that alternative can be used for testing if above breaks

  # For generating planetary system properties
  add_param_fixed(sim_param,"generate_planetary_system", generate_planetary_system_clustered)

  add_param_fixed(sim_param,"generate_num_clusters",generate_num_clusters_poisson) 
  add_param_fixed(sim_param,"generate_num_planets_in_cluster",generate_num_planets_in_cluster_poisson)
  add_param_active(sim_param,"log_rate_clusters",log(2.0))
  add_param_fixed(sim_param,"max_clusters_in_sys",10)
  add_param_active(sim_param,"log_rate_planets_per_cluster",log(3.0))
  add_param_fixed(sim_param,"max_planets_in_cluster",10)

  # generate_num_planets_in_cluster currently calls: generate_periods_power_law & generate_sizes_power_law
  add_param_fixed(sim_param,"power_law_P",0.3)
  add_param_fixed(sim_param,"power_law_r",-2.44)
  add_param_fixed(sim_param,"min_period",1.0)
  add_param_fixed(sim_param,"max_period",400.0)
  add_param_fixed(sim_param,"min_radius",0.5*ExoplanetsSysSim.earth_radius)
  add_param_fixed(sim_param,"max_radius",20.*ExoplanetsSysSim.earth_radius)

  # generate_num_planets_in_cluster currently use these for the Inclination distribution
  add_param_fixed(sim_param,"sigma_incl",3.0) # degrees; 0 = coplanar w/ generate_kepler_target_simple; ignored by generate_planetary_system_uncorrelated_incl
  add_param_fixed(sim_param,"sigma_incl_near_mmr",0.0)

  # generate_num_planets_in_cluster currently use these for the Eccentricity distribution
  add_param_fixed(sim_param,"generate_e_omega",ExoplanetsSysSim.generate_e_omega_rayleigh)
  add_param_fixed(sim_param,"sigma_hk",0.03)
  #add_param_fixed(sim_param,"sigma_hk_one",0.1)
  #add_param_fixed(sim_param,"sigma_hk_multi",0.03)

  # generate_num_planets_in_cluster currently use these for the Stability tests
  add_param_fixed(sim_param,"num_mutual_hill_radii",10.0)
  add_param_fixed(sim_param,"generate_planet_mass_from_radius",ExoplanetsSysSim.generate_planet_mass_from_radius_powerlaw)
  add_param_fixed(sim_param,"mr_power_index",2.0)
  add_param_fixed(sim_param,"mr_const",1.0)
  add_param_fixed(sim_param,"sigma_log_radius_in_cluster",1.0)
  add_param_fixed(sim_param,"sigma_logperiod_per_pl_in_cluster",1.0)

  # Functions to calculate observables from physical system properties
  add_param_fixed(sim_param,"calc_target_obs_single_obs",ExoplanetsSysSim.calc_target_obs_single_obs)   
  add_param_fixed(sim_param,"max_tranets_in_sys",8)     # SysSim ignores some planets in any systems with more than this many transiting planets to avoid wasting time on unphysical parameter values
  add_param_fixed(sim_param,"transit_noise_model",ExoplanetsSysSim.transit_noise_model_fixed_noise)
  #add_param_fixed(sim_param,"transit_noise_model",ExoplanetsSysSim.transit_noise_model_diagonal)
  # add_param_fixed(sim_param,"rng_seed",1234)   # If you want to be able to reproduce simulations
  add_param_fixed(sim_param,"read_target_obs",ExoplanetsSysSim.simulated_read_kepler_observations)  # Read Kepler observations to compare to from disk
  
  return sim_param
end


# Code for generating clustered planetary systems (once stable can move to src/planetary_system.jl)

function calc_hill_sphere(a::Float64, mu::Float64)
  a*(mu/3)^(1//3)
end

function calc_mutual_hill_radii{StarT<:StarAbstract}(ps::PlanetarySystem{StarT},pl1::Int64, pl2::Int64)
  mu = (ps.planet[pl1].mass+ps.planet[pl2].mass)/ps.star.mass
  a = 0.5*(ps.orbit[pl1].a+ps.orbit[pl2].a)
  calc_hill_sphere(a,mu)
end

function test_stability_circular(P::Vector{Float64},mass::Vector{Float64},star_mass::Float64, sim_param::SimParam)
   @assert length(P) == length(mass)
   const min_num_mutual_hill_radii = get_real(sim_param,"num_mutual_hill_radii")
   #println("# test_stability_circ: mass= ",mass," star_mass=",star_mass)
   #println("# test_stability_circ: mu= ",mass/star_mass, "; P= ",P)
   #print(  "# (aratio,Delta_H)= ")
   found_instability = false
   order = sortperm(P)
   a2 = semimajor_axis(P[order[1]],star_mass)
   for pl in 1:(length(P)-1)
       a1 = a2   # semimajor_axis(P[order[pl]],star_mass)
       a2 = semimajor_axis(P[order[pl+1]],star_mass)
       a = 0.5*(a1+a2)
       mu = (mass[order[pl]]+mass[order[pl+1]])/star_mass
       mutual_hill_radius = calc_hill_sphere(a,mu)
       #print("(",a2/a1,", ",(a2-a1)/mutual_hill_radius,"), ")
       if ! ( a2-a1  >= min_num_mutual_hill_radii*mutual_hill_radius )
          found_instability = true
          break
       end
   end # loop over neighboring planet pairs within cluster
   #println("# test_stability_circ = ",!found_instability)
   return !found_instability
end

function test_stability(P::Vector{Float64},mass::Vector{Float64},star_mass::Float64, sim_param::SimParam; ecc::Vector{Float64} = zeros(length(P)))
   @assert length(P) == length(mass) == length(ecc)
   const min_num_mutual_hill_radii = get_real(sim_param,"num_mutual_hill_radii")
   found_instability = false
   order = sortperm(P)
   a2 = semimajor_axis(P[order[1]],star_mass)
   for pl in 1:(length(P)-1)
       a1 = a2   # semimajor_axis(P[order[pl]],star_mass)
       a2 = semimajor_axis(P[order[pl+1]],star_mass)
       a = 0.5*(a1+a2)
       mu = (mass[order[pl]]+mass[order[pl+1]])/star_mass
       mutual_hill_radius = calc_hill_sphere(a,mu)
       e1 = ecc[order[pl]]
       e2 = ecc[order[pl+1]]
       if ! ( a2*(1-e2)-a1*(1+e1)  >= min_num_mutual_hill_radii*mutual_hill_radius )
          found_instability = true
          break
       end
   end # loop over neighboring planet pairs within cluster
   return !found_instability
end

function calc_if_near_resonance(P::Vector{Float64})
   @assert issorted(P)   # TODO: OPT: Could remove once know it is safe
   const resonance_width = 0.05   # TODO: FEATURE Make a model parameter?  Would need sim_param
   const resonance_width_factor = 1+resonance_width
   const period_ratios_to_check = [ 2.0, 1.5, 4/3, 5/4 ]
   result = falses(length(P))
   if length(P) >= 2
      for i in 1:(length(P)-1)
          for period_ratio in period_ratios_to_check
              if P[i]*period_ratio <= P[i+1] <= P[i]*period_ratio*resonance_width_factor
                 result[i] = true
                 result[i+1] = true
                 break
              end # if near resonance
          end # period_ratio
      end # planets
   end # at least two planets
   return result
end

function generate_planet_periods_sizes_masses_in_cluster( star::StarAbstract, sim_param::SimParam; n::Int64 = 1 )  # TODO: IMPORTANT: Make this function work and test before using for science
   @assert n>=1
   const generate_planet_mass_from_radius = get_function(sim_param,"generate_planet_mass_from_radius")

   if n==1
      R = ExoplanetsSysSim.generate_sizes_power_law(star,sim_param)[1]
      mass = map(r->generate_planet_mass_from_radius(r,sim_param),R)
      P = [1.0] # generate_periods_power_law(star,sim_param)
      return P, R, mass
   end

   # If reach here, then at least 2 planets in cluster

   mean_R = ExoplanetsSysSim.generate_sizes_power_law(star,sim_param)[1]
   const sigma_log_radius_in_cluster = get_real(sim_param,"sigma_log_radius_in_cluster")
   #println("# mean_R = ",mean_R," sigma_log_radius_in_cluster= ",sigma_log_radius_in_cluster)
   Rdist = LogNormal(log(mean_R),sigma_log_radius_in_cluster)
   R = rand(Rdist,n)

   #println("# Rp = ", R)
   mass = map(r->generate_planet_mass_from_radius(r,sim_param),R)
   #println("# mass = ", mass)

   const sigma_logperiod_per_pl_in_cluster = get_real(sim_param,"sigma_logperiod_per_pl_in_cluster")
   log_mean_P = 0.0 # log(generate_periods_power_law(star,sim_param))
   # Note: Currently, drawing all periods within a cluster at once and either keeping or rejecting the whole cluster
   #       Should we instead draw periods one at a time?
   Pdist = LogNormal(log_mean_P,sigma_logperiod_per_pl_in_cluster*n)
   local P
   found_good_periods = false
   attempts = 0
   max_attempts = 100  # Note: Should this be a parameter?
   while !found_good_periods && attempts<max_attempts
      attempts += 1
      P = rand(Pdist,n)
      if test_stability_circular(P,mass,star.mass,sim_param)
        found_good_periods = true
     end
   end # while trying to draw periods
   if !found_good_periods
      println("# Warning: Did not find a good set of periods, sizes and masses for one cluster.")
      return fill(NaN,n), R, mass  # Return NaNs for periods to indicate failed
   end
   return P, R, mass    # Note can also return earlier if only one planet in cluster or if fail to generate a good set of values
end

function generate_num_planets_in_cluster_poisson(s::Star, sim_param::SimParam)
  const lambda::Float64 = exp(get_real(sim_param,"log_rate_planets_per_cluster"))
  const max_planets_in_cluster::Int64 = get_int(sim_param,"max_planets_in_cluster")
  ExoplanetsSysSim.generate_num_planets_poisson(lambda,max_planets_in_cluster,min_planets=1)
end

function generate_num_clusters_poisson(s::Star, sim_param::SimParam)
  const lambda::Float64 = exp(get_real(sim_param,"log_rate_clusters"))
  const max_clusters_in_sys::Int64 = get_int(sim_param,"max_clusters_in_sys")
  ExoplanetsSysSim.generate_num_planets_poisson(lambda,max_clusters_in_sys)
end

# This version generates clustered planetary systems, adaptation of python code from Matthias He
function generate_planetary_system_clustered(star::StarAbstract, sim_param::SimParam; verbose::Bool = false)  # TODO: Make this function work and test before using for science
  # load functions to use for drawing parameters
  const generate_num_clusters = get_function(sim_param,"generate_num_clusters")
  const generate_num_planets_in_cluster = get_function(sim_param,"generate_num_planets_in_cluster")

  # Generate a set of periods, planet radii, and planet masses.
  attempt_system = 0
  max_attempts_system = 100
  local num_pl, Plist, Rlist, masslist
  valid_system = false
  while !valid_system && attempt_system <= max_attempts_system

     # First, generate number of clusters (to attempt) and planets (to attempt) in each cluster
     num_clusters = generate_num_clusters(star,sim_param)::Int64
     num_pl_in_cluster = map(x->generate_num_planets_in_cluster(star, sim_param)::Int64, 1:num_clusters)
     num_pl = sum(num_pl_in_cluster)

     if( num_pl==0 )
       return PlanetarySystem(star)
     end
     Plist = Array{Float64}(num_pl)
     Rlist = Array{Float64}(num_pl)
     masslist = Array{Float64}(num_pl)
     @assert num_pl_in_cluster[1] >= 1
     pl_start = 1
     pl_stop = 0
     for c in 1:num_clusters
         pl_stop  += num_pl_in_cluster[c]
         Plist_tmp, Rlist[pl_start:pl_stop], masslist[pl_start:pl_stop] = generate_planet_periods_sizes_masses_in_cluster(star, sim_param, n = num_pl_in_cluster[c])
         valid_cluster = !any(isnan.(Plist))
         valid_period_scale = false
         attempt_period_scale = 0
         max_attempts_period_scale = 100
         while !valid_period_scale && attempt_period_scale<=max_attempts_period_scale && valid_cluster
             attempt_period_scale += 1
             period_scale = ExoplanetsSysSim.generate_periods_power_law(star,sim_param)
             Plist[pl_start:pl_stop] = Plist_tmp .* period_scale
             if test_stability_circular(Plist[1:pl_stop],masslist[1:pl_stop],star.mass,sim_param)  # Note: Should we include eccentricities in this test?
                valid_period_scale = true
             end
         end  # while !valid_period_scale...
         if !valid_period_scale
            Plist[pl_start:pl_stop] = NaN
         end
         pl_start += num_pl_in_cluster[c]
     end # for c in 1:num_clusters
     if any(isnan.(Plist))  # If any loop failed to generate valid planets, it should set a NaN in the period list
        keep = .!(isnan.(Plist))  # Currently, keeping clusters that could be fit, rather than throwing out entirely and starting from scratch.  Is this a good idea?  Matthias tried the other approach in his python code.
        num_pl = sum(keep)
        Plist = Plist[keep]
        Rlist = Rlist[keep]
        masslist = masslist[keep]
     end
     valid_system = true
     #= Note: This would be for drawing each cluster separately and then accepting or rejecting the whole lot.
             By testing for stability before adding each cluster, this last test should be unnecessary.
     if test_stability(Plist,masslist,star.mass,sim_param)
        valid_system = true
     end
     =#
     attempt_system += 1
  end # while !valid_system...

  # Now assign orbits with given periods, sizes, and masses.
  const generate_e_omega = get_function(sim_param,"generate_e_omega")
  const sigma_ecc::Float64 = haskey(sim_param,"sigma_hk") ? get_real(sim_param,"sigma_hk") : 0.0
  const sigma_incl = deg2rad(get_real(sim_param,"sigma_incl"))
  const sigma_incl_near_mmr = deg2rad(get_real(sim_param,"sigma_incl_near_mmr"))
    pl = Array{Planet}(num_pl)
    orbit = Array{Orbit}(num_pl)
    incl_sys = acos(rand())
    idx = sortperm(Plist)       # TODO OPT: Check to see if sorting is significant time sink.  If so, could reduce redundant sortperm
    is_near_resonance = calc_if_near_resonance(Plist[idx])
    for i in 1:num_pl
      # if verbose   println("i=",i," idx=",idx," Plist=",Plist[idx] );     end
      if haskey(sim_param,"sigma_hk_one") && haskey(sim_param,"sigma_hk_multi")
         sigma_ecc = num_pl == 1 ? get_real(sim_param,"sigma_hk_one") : get_real(sim_param,"sigma_hk_multi")
      end
      (ecc,  omega) = generate_e_omega(sigma_ecc)::Tuple{Float64,Float64}  # WARNING: Not testing for stability after eccentricites drawn.
      sigma_incl_use = is_near_resonance[i] ? sigma_incl_near_mmr : sigma_incl
      incl_mut = sigma_incl_use*sqrt(randn()^2+randn()^2) # rand(Distributions.Rayleigh(sigma_incl))
      asc_node = 2pi*rand()
      mean_anom = 2pi*rand()
      incl =  incl_mut!=zero(incl_mut) ? acos( cos(incl_sys)*cos(incl_mut) + sin(incl_sys)*sin(incl_mut)*cos(asc_node) ) : incl_sys
      orbit[idx[i]] = Orbit(Plist[idx[i]],ecc,incl,omega,asc_node,mean_anom)
      pl[i] = Planet( Rlist[idx[i]],  masslist[idx[i]] )
    end # for i in 1:num_pl

  return PlanetarySystem(star,pl,orbit)
end

## summary_statistics borrowed from multiple_planets example, eventually should be merged into main branch

# Compile indices for N-tranet systems
function calc_summary_stats_idx_n_tranets!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  if haskey(css.cache,"idx_n_tranets")
     return css.cache["idx_n_tranets"]
  end
     max_tranets_in_sys = get_int(param,"max_tranets_in_sys")    
     idx_n_tranets = Vector{Int64}[ Int64[] for m = 1:max_tranets_in_sys]
     for n in 1:max_tranets_in_sys-1
       idx_n_tranets[n] = find(x::KeplerTargetObs-> length(x.obs) == n, cat_obs.target )
     end
     idx_n_tranets[max_tranets_in_sys] = find(x::KeplerTargetObs-> length(x.obs) >= max_tranets_in_sys, cat_obs.target )
     css.cache["idx_n_tranets"] = idx_n_tranets
  return idx_n_tranets 
end

# Count total number of tranets using lists of indices for N-tranet systems
function calc_summary_stats_num_tranets!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  if haskey(css.stat,"num_tranets")
     return css.stat["num_tranets"]
  elseif haskey(css.cache,"num_tranets")
     css.stat["num_tranets"] = css.cache["num_tranets"]
     return css.stat["num_tranets"]
  end
     idx_n_tranets = calc_summary_stats_idx_n_tranets!(css,cat_obs,param)
     num_tranets = 0
     for n in 1:length(idx_n_tranets)
         num_tranets += n*length(idx_n_tranets[n])
     end
     css.stat["num_tranets"] = num_tranets
  return num_tranets 
end

function calc_summary_stats_num_targets!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam ; trueobs_cat::Bool = false)
  if !trueobs_cat
    css.stat["num targets"] = get_int(param,"num_targets_sim_pass_one")
  else
    css.stat["num targets"] = get_int(param,"num_kepler_targets")
  end
end

function calc_summary_stats_num_n_tranet_systems!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  idx_n_tranets = calc_summary_stats_idx_n_tranets!(css,cat_obs,param)
  #max_tranets_in_sys = get_int(param,"max_tranets_in_sys")    
  num_n_tranet_systems = map(n->length(idx_n_tranets[n]), 1:length(idx_n_tranets) )
  #for n in 1:length(idx_n_tranets)
  #  num_n_tranet_systems[n] = length(idx_n_tranets[n])
  #end
  css.stat["num n-tranet systems"] = num_n_tranet_systems
end

function calc_summary_stats_duration_ratios_neighbors!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  if haskey(css.stat,"duration_ratio_list")
     return css.stat["duration_ratio_list"]
  elseif haskey(css.cache,"duration_ratio_list")
     return css.cache["duration_ratio_list"]
  end
  idx_n_tranets = calc_summary_stats_idx_n_tranets!(css,cat_obs,param)
  @assert length(idx_n_tranets) >= 1 

  # Calculate how many duration ratios there will be & allocate storage
  num_ratios = 0
  for i in 2:length(idx_n_tranets)
     num_ratios += length(idx_n_tranets[i])*(i-1)
  end
  duration_ratio_list = Array{Float64}(num_ratios)
 
  k = 0
  for n in 2:length(idx_n_tranets)         # Loop over number of tranets in system
    for i in idx_n_tranets[n]              # Loop over systems with n tranets
       period_in_sys = Array{Float64}(n)
       duration_in_sys = Array{Float64}(n)
       for j in 1:n                        # Loop over periods within a system
         period_in_sys[j] = cat_obs.target[i].obs[j].period
         duration_in_sys[j] = cat_obs.target[i].obs[j].duration
       end
       perm = sortperm(period_in_sys)
       for j in 1:(n-1)                       # Loop over period ratios within a system
          inv_period_ratio = period_in_sys[perm[j+1]]/period_in_sys[perm[j]]
          if 1<inv_period_ratio<Inf
             k += 1
             duration_ratio_list[k] = duration_in_sys[perm[j]]/duration_in_sys[perm[j+1]] * inv_period_ratio^(1//3)
          end
       end
    end
  end
  resize!(duration_ratio_list,k)
  css.stat["duration_ratio_list"] = duration_ratio_list

  return duration_ratio_list
end

function calc_summary_stats_period_radius_ratios_neighbors_internal!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  if haskey(css.stat,"period_ratio_list") && haskey(css.stat,"radius_ratio_list")
     return (css.stat["period_ratio_list"], css.stat["radius_ratio_list"])
  elseif haskey(css.cache,"period_ratio_list") && haskey(css.cache,"radius_ratio_list")
     return (css.cache["period_ratio_list"], css.cache["radius_ratio_list"])
  end
  idx_n_tranets = calc_summary_stats_idx_n_tranets!(css,cat_obs,param)
  @assert length(idx_n_tranets) >= 1 

  # Calculate how many period ratios there will be & allocate storage
  num_ratios = 0
  for i in 2:length(idx_n_tranets)
     num_ratios += length(idx_n_tranets[i])*(i-1)
  end
  period_ratio_list = Array{Float64}(num_ratios)
  radius_ratio_list = Array{Float64}(num_ratios)
 
  k = 0
  for n in 2:length(idx_n_tranets)         # Loop over number of tranets in system
    period_in_sys = Array{Float64}(n)
    #radius_in_sys = Array{Float64}(n)
    depth_in_sys = Array{Float64}(n)
    for i in idx_n_tranets[n]              # Loop over systems with n tranets
       for j in 1:n                        # Loop over periods within a system
         period_in_sys[j] = cat_obs.target[i].obs[j].period
         #radius_in_sys[j] = sqrt(cat_obs.target[i].obs[j].depth)*cat_obs.target[i].star.radius
         depth_in_sys[j] = cat_obs.target[i].obs[j].depth
       end
       perm = sortperm(period_in_sys)
       for j in 1:(n-1)                       # Loop over period ratios within a system
          k = k+1
          period_ratio_list[k] = period_in_sys[perm[j]]/period_in_sys[perm[j+1]]
          #radius_ratio_list[k] = radius_in_sys[perm[j]]/radius_in_sys[perm[j+1]]
          radius_ratio_list[k] = sqrt(depth_in_sys[perm[j]]/depth_in_sys[perm[j+1]]) 
       end
    end
  end
  css.cache["period_ratio_list"] = period_ratio_list
  css.cache["radius_ratio_list"] = radius_ratio_list

  return (period_ratio_list,radius_ratio_list)
end


function calc_summary_stats_period_radius_ratios_neighbors!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  (period_ratio_list,radius_ratio_list) = calc_summary_stats_period_radius_ratios_neighbors!(css,cat_obs,param)
  css.stat["period_ratio_list"] = period_ratio_list
  css.stat["radius_ratio_list"] = radius_ratio_list
  return (period_ratio_list, radius_ratio_list)
end
function calc_summary_stats_period_ratios_neighbors!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  period_ratio_list = calc_summary_stats_period_radius_ratios_neighbors_internal!(css,cat_obs,param)[1]
  css.stat["period_ratio_list"] = period_ratio_list
  return period_ratio_list
end
function calc_summary_stats_radius_ratios_neighbors!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  radius_ratio_list = calc_summary_stats_period_radius_ratios_neighbors_internal!(css,cat_obs,param)[2]
  css.stat["radius_ratio_list"] = radius_ratio_list
  return radius_ratio_list
end

function calc_summary_stats_mean_std_log_period_depth!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  # Allocate arrays to store values for each tranet
  num_tranets  = calc_summary_stats_num_tranets!(css, cat_obs, param)
  period_list = zeros(num_tranets)
  depth_list = zeros(num_tranets)
  #weight_list = ones(num_tranets)

  idx_n_tranets = calc_summary_stats_idx_n_tranets!(css, cat_obs, param)
  max_tranets_in_sys = get_int(param,"max_tranets_in_sys") 
  @assert max_tranets_in_sys >= 1
   i = 1   # tranet id
   for targ in cat_obs.target                        # For each target 
     for j in 1:min(length(targ.obs),max_tranets_in_sys)          # For each tranet around that target (but truncated if too many tranets in one system)
         #println("# i= ",i," j= ",j)
         period_list[i] = targ.obs[j].period
         depth_list[i] = targ.obs[j].depth
         #weight_list[i] = 1.0
         i = i+1
      end
   end

  css.cache["P list"] = period_list                                     # We can store whole lists, e.g., if we want to compute K-S distances
  css.cache["depth list"] = depth_list
  #css.cache["weight list"] = weight_list

  idx_good = Bool[ period_list[i]>0.0 && depth_list[i]>0.0 for i in 1:length(period_list) ]
  log_period_list = log10(period_list[idx_good])
  log_depth_list = log10(depth_list[idx_good])
  css.stat["mean log10 P"]  =  mean_log_P = mean(log_period_list)
  css.stat["mean log10 depth"]  =  mean_log_depth = mean(log_depth_list)
  css.stat["std log10 P"]  = std_log_P = stdm(log_period_list,mean_log_P)
  css.stat["std log10 depth"]  = std_log_depth = stdm(log_depth_list,mean_log_depth)

  return (mean_log_P, std_log_P, mean_log_depth, std_log_depth) 
end

function calc_summary_stats_cuml_period_depth!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  # Allocate arrays to store values for each tranet
  num_tranets  = calc_summary_stats_num_tranets!(css, cat_obs, param)
  period_list = zeros(num_tranets)
  depth_list = zeros(num_tranets)
  #weight_list = ones(num_tranets)

  idx_n_tranets = calc_summary_stats_idx_n_tranets!(css, cat_obs, param)
  max_tranets_in_sys = get_int(param,"max_tranets_in_sys") 
  @assert max_tranets_in_sys >= 1
  i = 0   # tranet id
   for targ in cat_obs.target                        # For each target 
     for j in 1:min(length(targ.obs),max_tranets_in_sys)          # For each tranet around that target (but truncated if too many tranets in one system)
         i = i+1
         #println("# i= ",i," j= ",j)
         period_list[i] = targ.obs[j].period
         depth_list[i] = targ.obs[j].depth
         #weight_list[i] = 1.0
      end
   end
  resize!(period_list,i)
  resize!(depth_list,i)
  css.stat["P list"] = period_list                                     # We can store whole lists, e.g., if we want to compute K-S distances
  css.stat["depth list"] = depth_list
  #css.cache["weight list"] = weight_list
  #println("# P list = ",period_list)
  return (period_list,depth_list)
end


function calc_summary_stats_obs_binned_rates!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam; trueobs_cat::Bool = false)
  num_tranets  = calc_summary_stats_num_tranets!(css, cat_obs, param)
  idx_n_tranets = calc_summary_stats_idx_n_tranets!(css, cat_obs, param)

  idx_tranets = find(x::KeplerTargetObs-> length(x.obs) > 0, cat_obs.target)::Array{Int64,1}             # Find indices of systems with at least 1 tranet = potentially detectable transiting planet
  css.cache["idx_tranets"] = idx_tranets                                   # We can save lists of indices to summary stats for pass 2, even though we will not use these for computing a distance or probability

  #=
  max_tranets_in_sys = get_int(param,"max_tranets_in_sys")    # Demo that simulation parameters can specify how to evalute models, too
  @assert max_tranets_in_sys >= 1

  # Count total number of tranets and compile indices for N-tranet systems
  num_tranets = 0
  idx_n_tranets = Vector{Int64}[ Int64[] for m = 1:max_tranets_in_sys]
  for n in 1:max_tranets_in_sys-1
    idx_n_tranets[n] = find(x::KeplerTargetObs-> length(x.obs) == n, cat_obs.target )
    num_tranets += n*length(idx_n_tranets[n])
  end
  idx_n_tranets[max_tranets_in_sys] = find(x::KeplerTargetObs-> length(x.obs) >= max_tranets_in_sys, cat_obs.target )
  css.cache["idx_n_tranets"] = idx_n_tranets

  num_tranets += max_tranets_in_sys*length(idx_n_tranets[max_tranets_in_sys])  # WARNING: this means we need to ignore planets w/ indices > max_tranets_in_sys
  if ( length( find(x::KeplerTargetObs-> length(x.obs) > max_tranets_in_sys, cat_obs.target ) ) > 0)   # Make sure max_tranets_in_sys is at least big enough for observed systems
    warn("Observational data has more transiting planets in one systems than max_tranets_in_sys allows.")
  end
  num_tranets  = convert(Int64,num_tranets)            # TODO OPT: Figure out why is not this already an Int.  I may be doing something that prevents some optimizations
  css.cache["num_tranets"] = num_tranets                                   

  num_sys_tranets = zeros(max_tranets_in_sys)                           # Since observed data, don't need to calculate probabilities.
  for n in 1:max_tranets_in_sys                                         # Make histogram of N-tranet systems
    num_sys_tranets[n] = length(idx_n_tranets[n])
  end
  css.stat["num_sys_tranets"] = num_sys_tranets
  css.stat["planets detected"] = num_tranets 
  =#

  period_list = zeros(num_tranets)
  radius_list = zeros(num_tranets)
  weight_list = ones(num_tranets)

  n = 1    # tranet id
  for i in 1:length(cat_obs.target)
    for j in 1:num_planets(cat_obs.target[i])
      period_list[n] = cat_obs.target[i].obs[j].period
      radius_list[n] = sqrt(cat_obs.target[i].obs[j].depth)*cat_obs.target[i].star.radius
      n = n+1
    end
  end

  limitP::Array{Float64,1} = get_any(param, "p_lim_arr", Array{Float64,1})
  limitRp::Array{Float64,1} = get_any(param, "r_lim_arr", Array{Float64,1})
  @assert length(limitP)>=2 && length(limitRp)>=2

  np_bin = zeros((length(limitP)-1) * (length(limitRp)-1))
  np_bin_idx = 1
  for i in 1:(length(limitP)-1)
    P_match = find(x -> ((x > limitP[i]) && (x < limitP[i+1])), period_list)
    for j in 1:(length(limitRp)-1)
      R_match = find(x -> ((x > limitRp[j]) && (x < limitRp[j+1])), radius_list)
      
      bin_match = intersect(P_match, R_match)

      np_bin[np_bin_idx] = sum(weight_list[bin_match])
      np_bin_idx += 1
    end
  end

  #css.stat["planets detected"] = sum(np_bin)
  css.stat["planets table"] = np_bin

  return css
end


function calc_summary_stats_model(cat_obs::KeplerObsCatalog, param::SimParam; trueobs_cat::Bool = false)
  css = CatalogSummaryStatistics()
  calc_summary_stats_num_targets!(css,cat_obs,param,trueobs_cat=trueobs_cat)
  calc_summary_stats_num_tranets!(css,cat_obs,param)
  calc_summary_stats_num_n_tranet_systems!(css,cat_obs,param)
  calc_summary_stats_cuml_period_depth!(css,cat_obs,param)
  #calc_summary_stats_obs_binned_rates!(css,cat_obs,param)
  #calc_summary_stats_mean_std_log_period_depth!(css,cat_obs,param)
  calc_summary_stats_period_ratios_neighbors!(css,cat_obs,param)
  calc_summary_stats_duration_ratios_neighbors!(css,cat_obs,param)
  return css
end

## abc_distance
function calc_distance_model_ks(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  d1 = calc_distance_num_planets(summary1,summary2,sim_param,verbose=verbose)
  d2 = calc_distance_num_n_tranet_systems(summary1,summary2,sim_param,verbose=verbose)
  #d3 = calc_distance_num_planets_binned(summary1,summary2,sim_param,verbose=verbose)
  #d4 = calc_distance_mean_std_log_period_depth(summary1,summary2,sim_param,verbose=verbose)
  d5 = calc_distance_ks_period(summary1,summary2,sim_param,verbose=verbose)
  d6 = calc_distance_ks_depth(summary1,summary2,sim_param,verbose=verbose)
  d7 = calc_distance_ks_period_ratios(summary1,summary2,sim_param,verbose=verbose)
  d8 = calc_distance_ks_duration_ratios(summary1,summary2,sim_param,verbose=verbose)
  return vcat(d1, d2, d5, d6, d7, d8)
end

function calc_distance_model_kl(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  d1 = calc_distance_kl_num_planets(summary1,summary2,sim_param,verbose=verbose)
  d2 = calc_distance_kl_num_n_tranet_systems(summary1,summary2,sim_param,verbose=verbose)
  #d2 = calc_distance_hellinger_num_n_tranet_systems(summary1,summary2,sim_param,verbose=verbose)
  #d3 = calc_distance_num_planets_binned(summary1,summary2,sim_param,verbose=verbose)
  #d4 = calc_distance_mean_std_log_period_depth(summary1,summary2,sim_param,verbose=verbose)
  d5 = calc_distance_kl_period(summary1,summary2,sim_param,verbose=verbose)
  d6 = calc_distance_kl_depth(summary1,summary2,sim_param,verbose=verbose)
  d7 = calc_distance_kl_period_ratios(summary1,summary2,sim_param,verbose=verbose)
  d8 = calc_distance_kl_duration_ratios(summary1,summary2,sim_param,verbose=verbose)
  return vcat(d1, d2, d5, d6, d7, d8)
end
calc_distance_model(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false) = calc_distance_model_kl(summary1,summary2,sim_param,verbose=verbose)

function calc_distance_num_planets(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
    np1 = haskey(summary1.stat,"num_tranets") ? summary1.stat["num_tranets"] : summary1.stat["expected planets detected"]
    np2 = haskey(summary2.stat,"num_tranets") ? summary2.stat["num_tranets"] : summary2.stat["expected planets detected"]
    #println("np1 = ",np1,", np2 = ",np2)

    dist_np = dist_L1_abs(np1/summary1.stat["num targets"], np2/summary2.stat["num targets"])
    #println("np1 (normalized) = ",np1/summary1.stat["num targets"],", np2 (normalized) = ",np2/summary2.stat["num targets"],", d = ",dist_np)
  return dist_np
end

function calc_distance_num_planets_binned(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
    np1 = haskey(summary1.stat,"planets table") ? summary1.stat["planets table"] : summary1.stat["expected planets table"]
    np2 = haskey(summary2.stat,"planets table") ? summary2.stat["planets table"] : summary2.stat["expected planets table"]
    #println("np1 = ",np1,", np2 = ",np2)

    dist_np_bin = zeros(length(np1))
    for n in 1:length(np1)
      dist_np_bin[n] = dist_L1_abs(np1[n]/summary1.stat["num targets"], np2[n]/summary2.stat["num targets"])
      #println("True # [Bin ", n,"] = ",np1[n],", Expected # [Bin ", n,"] = ",np2[n])
    end
    #println("np1 (normalized) = ",np1/summary1.stat["num targets"],", np2 (normalized) = ",np2/summary2.stat["num targets"],", dist = ",dist_np_bin)
  return dist_np_bin
end

function calc_distance_num_n_tranet_systems(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  max_tranets_in_sys = get_int(sim_param,"max_tranets_in_sys")    
  d = zeros(max_tranets_in_sys)    
  for n in 1:max_tranets_in_sys   
    d[n] = dist_L1_abs(summary1.stat["num n-tranet systems"][n]/summary1.stat["num targets"], summary2.stat["num n-tranet systems"][n]/summary2.stat["num targets"])
  end
  return d
end

function calc_distance_mean_std_log_period_depth(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
    d1 = dist_L1_abs(summary1.stat["mean log10 P"], summary2.stat["mean log10 P"])
    d2 = dist_L1_abs(summary1.stat["mean log10 depth"], summary2.stat["mean log10 depth"])
    d3 = dist_L1_abs(summary1.stat["std log10 P"], summary2.stat["std log10 P"])
    d4 = dist_L1_abs(summary1.stat["std log10 depth"], summary2.stat["std log10 depth"])
    return [d1, d2, d3, d4]
end

# compute supremum of differences between empirical cdfs.
# Borrowed from JuliaStats/HypothesisTests.jl
function ksstats{T<:Real, S<:Real}(x::AbstractVector{T}, y::AbstractVector{S})
    n_x, n_y = length(x), length(y)
    sort_idx = sortperm([x; y])
    pdf_diffs = [ones(n_x)/n_x; -ones(n_y)/n_y][sort_idx]
    cdf_diffs = cumsum(pdf_diffs)
    deltap = maximum(cdf_diffs)
    deltan = -minimum(cdf_diffs)
    delta = max(deltap,deltan)
    (n_x, n_y, deltap, deltan, delta)
end

function calc_distance_ks_period(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  samp1 = summary1.stat["P list"] 
  samp2 = summary2.stat["P list"] 
  return ksstats(samp1,samp2)[5]
end

function calc_distance_ks_depth(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  samp1 = summary1.stat["depth list"] 
  samp2 = summary2.stat["depth list"] 
  return ksstats(samp1,samp2)[5]
end

function calc_distance_ks_period_ratios(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  samp1 = summary1.stat["period_ratio_list"] 
  samp2 = summary2.stat["period_ratio_list"] 
  return ksstats(samp1,samp2)[5]
end

function calc_distance_ks_duration_ratios(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  samp1 = summary1.stat["duration_ratio_list"] 
  samp2 = summary2.stat["duration_ratio_list"] 
  return ksstats(samp1,samp2)[5]
end

# Function for Relative Entropy / K-L divergence
include("kde.jl")

function calc_distance_kl_num_planets(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
    np1 = (haskey(summary1.stat,"num_tranets") ? summary1.stat["num_tranets"] : summary1.stat["expected planets detected"]) 
    np2 = (haskey(summary2.stat,"num_tranets") ? summary2.stat["num_tranets"] : summary2.stat["expected planets detected"]) 
    ntarg1 = summary1.stat["num targets"]
    ntarg2 = summary2.stat["num targets"]
    mintarg = min(ntarg1,ntarg2)
    alpha1 = 1+np1
    beta1 = 1+ntarg1
    #mu1 = alpha1/beta1
    #sigma1 = mu1 /sqrt(alpha1)
    #dist = alpha+lgamma(alpha)+(1-alpha)*digamma(alpha) -log(beta)
    alpha2 = 1+np2
    beta2 = 1+ntarg2
    #dist -= alpha+lgamma(alpha)+(1-alpha)*digamma(alpha) -log(beta)
    #mu2 = alpha2/beta2
    #sigma2 = mu2 /sqrt(alpha1)
    #dist = 0.5*( (sigma2/sigma1)^2 - 1 + 2*log(sigma1/sigma2) + (mu2-mu1)^2/sigma1^2 )  # Normal approximation, substitute denom for difference in rates
    #= Gamma PDF
    dist = abs( (alpha2-alpha1)*digamma(alpha1) - lgamma(alpha1) + lgamma(alpha2)  ) # Gamma approximation, same beta, why need abs? 
    dist += abs( (alpha1-alpha2)*digamma(alpha2) - lgamma(alpha2) + lgamma(alpha1)  )  # Gamma approximation, same beta, why need abs? 
    print("# a1=",alpha1, " a2=",alpha2, " digamma=",digamma(alpha1), " ", digamma(alpha2)," gamma term=",lgamma(alpha1)-lgamma(alpha2))
    if ntarg1!=ntarg2
       print("# ntarg=",ntarg1," ",ntarg2," dist(before betas)=", dist) 
       dist += alpha2*(log(beta1)-log(beta2)) + alpha1*(beta2-beta1)/beta1  # Gammma approximation, beta terms
       dist += alpha1*(log(beta2)-log(beta1)) + alpha2*(beta1-beta2)/beta2  # Gammma approximation, beta terms
       println(" dist (after betas)=",dist) 
    end
    =#
    ## #= Poisson 
    #dist = abs( alpha2-alpha1+alpha1*log(alpha1/alpha2) )
    #dist += abs( alpha2-alpha1+alpha1*log(alpha1/alpha2) )
    rate1 = alpha1/beta1*mintarg
    rate2 = alpha2/beta2*mintarg
    dist = (rate1-rate2)*log(rate1/rate2) 
    ## =#
    dist /= 2
    return dist
end

function calc_distance_kl_num_n_tranet_systems(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  max_tranets_in_sys = get_int(sim_param,"max_tranets_in_sys")   
  #= categorical distribution
  f1sum = sum(summary1.stat["num n-tranet systems"]) # summary1.stat["num targets"]
  f2sum = sum(summary2.stat["num n-tranet systems"]) # summary2.stat["num targets"]
  if !(f1sum>0 && f2sum>0) return 0.0  end 
  d = zeros(max_tranets_in_sys)
  for n in 1:max_tranets_in_sys   
    f1 = summary1.stat["num n-tranet systems"][n]/f1sum
    f2 = summary2.stat["num n-tranet systems"][n]/f2sum
    m = (f1+f2)/2
    if m>zero(m)
    if f1>zero(f1) 
       d[n] += 0.5*f1*log(f1/m)
    end
    if f2>zero(f2) 
       d[n] += 0.5*f2*log(f2/m)
    end
#   else
#       d += Inf
    end
  end
  =#
  # Poisson distributions for each
    ntarg1 = summary1.stat["num targets"]
    ntarg2 = summary2.stat["num targets"]
    #mintarg = min(ntarg1,ntarg2)
  d = zeros(max_tranets_in_sys)
  for n in 1:max_tranets_in_sys   
    f1 = (1+summary1.stat["num n-tranet systems"][n])/(1+ntarg1) # *(1+mintarg)
    f2 = (1+summary2.stat["num n-tranet systems"][n])/(1+ntarg2) # *(1+mintarg)
    d[n] = (f1-f2)*log(f1/f2)
  end
  return d
end
function calc_distance_hellinger_num_n_tranet_systems(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  max_tranets_in_sys = get_int(sim_param,"max_tranets_in_sys")    
  f1sum = sum(summary1.stat["num n-tranet systems"]) # summary1.stat["num targets"]
  f2sum = sum(summary2.stat["num n-tranet systems"]) # summary2.stat["num targets"]
  d = 1
  for n in 1:max_tranets_in_sys   
    f1 = summary1.stat["num n-tranet systems"][n]/f1sum
    f2 = summary2.stat["num n-tranet systems"][n]/f2sum
    d -= sqrt(f1*f2)
  end
  #d = sqrt(d)
  return d
end

function calc_distance_kl_period(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  #println("# P list 1: n=",length(summary1.stat["P list"])," min=",minimum(summary1.stat["P list"]), " max=",maximum(summary1.stat["P list"]))
  #println("# P list 2: n=",length(summary2.stat["P list"])," min=",minimum(summary2.stat["P list"]), " max=",maximum(summary2.stat["P list"]))
  samp1 = log(summary1.stat["P list"] )
  samp2 = log(summary2.stat["P list"] )
  calc_kl_distance_ab(samp1,samp2,log(0.5),log(320.) )
end

function calc_distance_kl_depth(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  samp1 = log(summary1.stat["depth list"] )
  samp2 = log(summary2.stat["depth list"] )
  calc_kl_distance_ab(samp1,samp2,log(0.000025),log(0.025) )
end

function calc_distance_kl_period_ratios(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  min_ratios_to_compute_distances = 3
  distance_when_not_enough_ratios = 10.0
  samp1 = summary1.stat["period_ratio_list"] 
  samp2 = summary2.stat["period_ratio_list"] 
  if length(samp1)<min_ratios_to_compute_distances || length(samp2)<min_ratios_to_compute_distances
     return distance_when_not_enough_ratios
  end
  calc_kl_distance_ab(samp1,samp2,0.0,1.0)
end

function calc_distance_kl_duration_ratios(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  min_ratios_to_compute_distances = 3
  distance_when_not_enough_ratios = 10.0
  samp1 = log( summary1.stat["duration_ratio_list"] )
  samp2 = log( summary2.stat["duration_ratio_list"] )
  if length(samp1)<min_ratios_to_compute_distances || length(samp2)<min_ratios_to_compute_distances
     return distance_when_not_enough_ratios
  end
  calc_kl_distance_ab(samp1,samp2,-log(10.0),log(10.0) )
end


## Old code for generating stellar properties # TODO: WARNING: Should eventually use version in main branch to make sure have recent improvements

## stellar_table
function setup_star_table_christiansen(sim_param::SimParam; force_reread::Bool = false)
  #global df
  df = ExoplanetsSysSim.StellarTable.df
  if haskey(sim_param,"read_stellar_catalog") && !force_reread
     return df
     #return data
  end
  stellar_catalog_filename = convert(ASCIIString,joinpath(Pkg.dir("ExoplanetsSysSim"), "data", convert(ASCIIString,get(sim_param,"stellar_catalog","q1_q17_dr25_stellar.csv")) ) )
  df = setup_star_table_christiansen(stellar_catalog_filename)
  add_param_fixed(sim_param,"read_stellar_catalog",true)
  set_star_table(df)
  return df  
end

function setup_star_table_christiansen(filename::ASCIIString; force_reread::Bool = false)
  #global df, usable
  df = ExoplanetsSysSim.StellarTable.df
  usable = ExoplanetsSysSim.StellarTable.usable
  if ismatch(r".jld$",filename)
  try 
    data = load(filename)
    df::DataFrame = data["stellar_catalog"]
    usable::Array{Int64,1} = data["stellar_catalog_usable"]
    set_star_table(df, usable)
  catch
    error(string("# Failed to read stellar catalog >",filename,"< in jld format."))
  end
  else
  try 
    df = readtable(filename)
  catch
    error(string("# Failed to read stellar catalog >",filename,"< in ascii format."))
  end

  has_mass = ! (isna(df[:mass]) | isna(df[:mass_err1]) | isna(df[:mass_err2]))
  has_radius = ! (isna(df[:radius]) | isna(df[:radius_err1]) | isna(df[:radius_err2]))
  has_dens = ! (isna(df[:dens]) | isna(df[:dens_err1]) | isna(df[:dens_err2]))
  has_rest = ! (isna(df[:rrmscdpp04p5]) | isna(df[:dataspan]) | isna(df[:dutycycle]))
  in_Q1Q12 = []
  for x in df[:st_quarters]
    subx = string(x)
    subx = ("0"^(17-length(subx)))*subx
    indQ = search(subx, '1')
    if ((indQ < 1) | (indQ > 12))
      push!(in_Q1Q12, false)
    else
      push!(in_Q1Q12, true)
    end
  end
  is_FGK = []
  for x in 1:length(df[:teff])
    if ((df[x,:teff] > 4000.0) & (df[x,:teff] < 7000.0) & (df[x,:logg] > 4.0))
      push!(is_FGK, true)
    else
      push!(is_FGK, false)
    end
  end
  is_usable = has_radius & is_FGK & has_mass & has_rest #& in_Q1Q12 # & has_dens
  #if contains(filename,"q1_q12_christiansen.jld")
  if contains(filename,"q1_q12_christiansen")   # TODO: Ask Danely what he's trying to do here.
    is_usable = is_usable & in_Q1Q12
  end
  # See options at: http://exoplanetarchive.ipac.caltech.edu/docs/API_keplerstellar_columns.html
  # TODO SCI DETAIL or IMPORTANT?: Read in all CDPP's, so can interpolate?
  symbols_to_keep = [ :kepid, :mass, :mass_err1, :mass_err2, :radius, :radius_err1, :radius_err2, :dens, :dens_err1, :dens_err2, :rrmscdpp04p5, :dataspan, :dutycycle ]
  delete!(df, [~(x in symbols_to_keep) for x in names(df)])    # delete columns that we won't be using anyway
  usable = find(is_usable)
  df = df[usable, symbols_to_keep]
  set_star_table(df, usable)
  end
  return df
end


## eval_model   # TODO: Update to test code in this model
function test_model()
  global sim_param_closure = setup_sim_param_model()
  cat_phys = generate_kepler_physical_catalog(sim_param_closure)
  cat_obs = observe_kepler_targets_single_obs(cat_phys,sim_param_closure)
  global summary_stat_ref_closure = calc_summary_stats_model(cat_obs,sim_param_closure)
  #global cat_phys_try_closure = generate_kepler_physical_catalog(sim_param_closure)
  #global cat_obs_try_closure  = observe_kepler_targets_single_obs(cat_phys_try_closure,sim_param_closure)
  #global summary_stat_try_closure  = calc_summary_stats_model(cat_obs_try_closure,cat_phys_try_closure,sim_param_closure)
  #param_guess = make_vector_of_sim_param(sim_param_closure)
  #evaluate_model_scalar_ret( param_guess)
end


