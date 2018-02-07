if !isdefined(:ExoplanetsSysSim) using ExoplanetsSysSim end 

## simulation_parameters
function setup_sim_param_model(args::Vector{String} = Array{String}(0) )   # allow this to take a list of parameter (e.g., from command line)
  sim_param = SimParam()
  # How many tatrges to generate
  #add_param_fixed(sim_param,"num_targets_sim_pass_one",150969)                      # Note this is used for the number of stars in the simulations, not necessarily related to number of Kepler targets
  add_param_fixed(sim_param,"num_targets_sim_pass_one",150960)                      # Note this is used for the number of stars in the simulations, not necessarily related to number of Kepler targets
  add_param_fixed(sim_param,"num_kepler_targets",150969)                            # Note this is used for the number of Kepler targets for the observational catalog

  # For generating target star properties
  add_param_fixed(sim_param,"generate_kepler_target",generate_kepler_target_from_table)
  add_param_fixed(sim_param,"star_table_setup",setup_star_table_christiansen)
#add_param_fixed(sim_param,"stellar_catalog","q1_q17_dr24_stellar.jld")
  add_param_fixed(sim_param,"stellar_catalog","q1_q17_dr24_stellar.csv")
  # add_param_fixed(sim_param,"generate_kepler_target",ExoplanetsSysSim.generate_kepler_target_simple)  # An alternative that alternative can be used for testing if above breaks

  # For generating planetary system properties
  add_param_fixed(sim_param,"generate_planetary_system", generate_planetary_system_clustered)

  add_param_fixed(sim_param,"generate_num_clusters",generate_num_clusters_poisson) 
  add_param_fixed(sim_param,"generate_num_planets_in_cluster",generate_num_planets_in_cluster_poisson)
  add_param_active(sim_param,"log_rate_clusters",log(0.45))
  add_param_fixed(sim_param,"max_clusters_in_sys",5)
  add_param_active(sim_param,"log_rate_planets_per_cluster",log(1.8))
  add_param_fixed(sim_param,"max_planets_in_cluster",10)

  # generate_num_planets_in_cluster currently calls: generate_periods_power_law & generate_sizes_power_law
  add_param_fixed(sim_param,"power_law_P",0.07)
  add_param_fixed(sim_param,"power_law_r",-1.54)
  add_param_fixed(sim_param,"min_period",5.0)
  add_param_fixed(sim_param,"max_period",300.0)
  add_param_fixed(sim_param,"min_radius",0.5*ExoplanetsSysSim.earth_radius)
  add_param_fixed(sim_param,"max_radius",10.*ExoplanetsSysSim.earth_radius)

  # generate_num_planets_in_cluster currently use these for the Inclination distribution
  add_param_fixed(sim_param,"sigma_incl",1.52) # degrees; 0 = coplanar w/ generate_kepler_target_simple; ignored by generate_planetary_system_uncorrelated_incl
  add_param_fixed(sim_param,"sigma_incl_near_mmr",0.0)

  # generate_num_planets_in_cluster currently use these for the Eccentricity distribution
  add_param_fixed(sim_param,"generate_e_omega",ExoplanetsSysSim.generate_e_omega_rayleigh)
  add_param_fixed(sim_param,"sigma_hk",0.05)
  #add_param_fixed(sim_param,"sigma_hk_one",0.1)
  #add_param_fixed(sim_param,"sigma_hk_multi",0.03)

  # generate_num_planets_in_cluster currently use these for the Stability tests
  add_param_fixed(sim_param,"num_mutual_hill_radii",15.0)
  add_param_fixed(sim_param,"generate_planet_mass_from_radius",ExoplanetsSysSim.generate_planet_mass_from_radius_powerlaw)
  add_param_fixed(sim_param,"mr_power_index",2.0)
  add_param_fixed(sim_param,"mr_const",1.0)
  add_param_fixed(sim_param,"mr_max_mass",1e3*ExoplanetsSysSim.earth_mass)
  add_param_fixed(sim_param,"sigma_log_radius_in_cluster",1.0)
  add_param_fixed(sim_param,"sigma_logperiod_per_pl_in_cluster",0.21)

  # Functions to calculate observables from physical system properties
  add_param_fixed(sim_param,"calc_target_obs_single_obs",ExoplanetsSysSim.calc_target_obs_single_obs)   
  add_param_fixed(sim_param,"max_tranets_in_sys",8)     # SysSim ignores some planets in any systems with more than this many transiting planets to avoid wasting time on unphysical parameter values
  add_param_fixed(sim_param,"transit_noise_model",ExoplanetsSysSim.transit_noise_model_fixed_noise)
  #add_param_fixed(sim_param,"transit_noise_model",ExoplanetsSysSim.transit_noise_model_diagonal)
  # add_param_fixed(sim_param,"rng_seed",1234)   # If you want to be able to reproduce simulations
  add_param_fixed(sim_param,"read_target_obs",ExoplanetsSysSim.simulated_read_kepler_observations)  # Read Kepler observations to compare to from disk
  
  # Read any customizations, so its easy to keep them separate from the essentials
  #=if isfile("param_custom.jl")
     include("param_custom.jl")
  end
  =#
  if isdefined(:add_param_custom)
     sim_param = add_param_custom(sim_param)
  end
  return sim_param
end

function test_setup_sim_param()
  setup_sim_param_model()
end


