if !isdefined(:ExoplanetsSysSim) using ExoplanetsSysSim end 
import Compat: UTF8String, ASCIIString

## simulation_parameters
function setup_sim_param_model(args::Vector{String} = Array{String}(0) )   # allow this to take a list of parameter (e.g., from command line)
  sim_param = SimParam()
  # How many tatrges to generate
  #add_param_fixed(sim_param,"num_targets_sim_pass_one",150061)                      # Note this is used for the number of stars in the simulations, not necessarily related to number of Kepler targets
  add_param_fixed(sim_param,"num_targets_sim_pass_one",78005)                      # Note this is used for the number of stars in the simulations, not necessarily related to number of Kepler targets
  add_param_fixed(sim_param,"num_kepler_targets",78005)                            # Note this is used for the number of Kepler targets for the observational catalog

  # For generating target star properties
  add_param_fixed(sim_param,"generate_kepler_target",generate_kepler_target_from_table)
  add_param_fixed(sim_param,"star_table_setup",setup_star_table_christiansen)
  #add_param_fixed(sim_param,"stellar_catalog","q1_q17_dr25_stellar.jld") #currently the JLD files do not have all the necessary fields
  add_param_fixed(sim_param,"stellar_catalog","q1q17_dr25_gaia_fgk.jld") #"q1_q17_dr25_stellar.csv"
  #add_param_fixed(sim_param,"generate_kepler_target",ExoplanetsSysSim.generate_kepler_target_simple)  # An alternative that alternative can be used for testing if above breaks

  # For generating planetary system properties
  add_param_fixed(sim_param,"generate_planetary_system", generate_planetary_system_clustered)

  add_param_fixed(sim_param,"generate_num_clusters",generate_num_clusters_poisson) 
  add_param_fixed(sim_param,"generate_num_planets_in_cluster",generate_num_planets_in_cluster_poisson)
  add_param_active(sim_param,"log_rate_clusters",log(3.0))
  add_param_fixed(sim_param,"max_clusters_in_sys",10)
  add_param_active(sim_param,"log_rate_planets_per_cluster",log(3.0))
  add_param_fixed(sim_param,"max_planets_in_cluster",10)

  # generate_num_planets_in_cluster currently calls: generate_periods_power_law
  add_param_fixed(sim_param,"generate_sizes",ExoplanetsSysSim.generate_sizes_broken_power_law) # To choose the way we draw planetary radii; if "generate_sizes_power_law", then takes "power_law_r"; if "generate_sizes_broken_power_law", then takes "power_law_r1", "power_law_r2", and "break_radius"
  add_param_active(sim_param,"power_law_P",0.5)
  add_param_fixed(sim_param,"power_law_r",-2.5)
  add_param_active(sim_param,"power_law_r1",-2.0)
  add_param_active(sim_param,"power_law_r2",-4.0)
  add_param_fixed(sim_param,"min_period",3.0)
  add_param_fixed(sim_param,"max_period",300.0)
  add_param_fixed(sim_param,"min_radius",0.5*ExoplanetsSysSim.earth_radius)
  add_param_fixed(sim_param,"max_radius",10.*ExoplanetsSysSim.earth_radius)
  add_param_active(sim_param,"break_radius",3.0*ExoplanetsSysSim.earth_radius)

  # generate_num_planets_in_cluster currently use these for the Inclination distribution
  add_param_active(sim_param,"sigma_incl",1.5) # degrees; 0 = coplanar w/ generate_kepler_target_simple; ignored by generate_planetary_system_uncorrelated_incl
  add_param_active(sim_param,"sigma_incl_near_mmr",1.5)

  add_param_fixed(sim_param,"max_incl_sys",0.0) #degrees; gives system inclinations from "max_incl_sys" (deg) to 90 (deg), so set to 0 for isotropic distribution of system inclinations; NOTE: make sure the difference between this and 90 (deg) is at least greater than "sigma_incl" and "sigma_incl_near_mmr"!

  # generate_num_planets_in_cluster currently use these for the Eccentricity distribution
  add_param_fixed(sim_param,"generate_e_omega",ExoplanetsSysSim.generate_e_omega_rayleigh)
  add_param_active(sim_param,"sigma_hk",0.05)
  #add_param_fixed(sim_param,"sigma_hk_one",0.1)
  #add_param_fixed(sim_param,"sigma_hk_multi",0.03)

  # generate_num_planets_in_cluster currently use these for the Stability tests
  add_param_active(sim_param,"num_mutual_hill_radii",8.0) #10.0
  add_param_fixed(sim_param,"generate_planet_mass_from_radius",generate_planet_mass_from_radius_Ning2018) # "ExoplanetsSysSim.generate_planet_mass_from_radius_powerlaw" or "generate_planet_mass_from_radius_Ning2018"
  add_param_fixed(sim_param,"mr_power_index",2.0)
  add_param_fixed(sim_param,"mr_const",1.0)
  add_param_fixed(sim_param,"mr_max_mass",1e3*ExoplanetsSysSim.earth_mass)
  add_param_active(sim_param,"sigma_log_radius_in_cluster",0.25)
  add_param_active(sim_param,"sigma_logperiod_per_pl_in_cluster",0.15)

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



function write_model_params(f, sim_param::SimParam)
    #This function writes all the model parameters to a file f as a header
    println(f, "# num_targets_sim_pass_one: ", get_int(sim_param,"num_targets_sim_pass_one"))
    println(f, "# max_incl_sys: ", get_real(sim_param,"max_incl_sys"))
    println(f, "# log_rate_clusters: ", get_real(sim_param,"log_rate_clusters"))
    println(f, "# max_clusters_in_sys: ", get_int(sim_param,"max_clusters_in_sys"))
    println(f, "# log_rate_planets_per_cluster: ", get_real(sim_param,"log_rate_planets_per_cluster"))
    println(f, "# max_planets_in_clusters: ", get_int(sim_param,"max_planets_in_cluster"))
    println(f, "# power_law_P: ", get_real(sim_param,"power_law_P"))
    println(f, "# min_period: ", get_real(sim_param,"min_period"))
    println(f, "# max_period: ", get_real(sim_param,"max_period"))

    if string(get_function(sim_param,"generate_sizes")) == "ExoplanetsSysSim.generate_sizes_power_law"
        println(f, "# power_law_r: ", get_real(sim_param,"power_law_r"))
    elseif string(get_function(sim_param,"generate_sizes")) == "ExoplanetsSysSim.generate_sizes_broken_power_law"
        println(f, "# power_law_r1: ", get_real(sim_param,"power_law_r1"))
        println(f, "# power_law_r2: ", get_real(sim_param,"power_law_r2"))
        println(f, "# break_radius (R_earth): ", get_real(sim_param,"break_radius")/ExoplanetsSysSim.earth_radius)
    end

    println(f, "# min_radius (R_earth): ", get_real(sim_param,"min_radius")/ExoplanetsSysSim.earth_radius)
    println(f, "# max_radius (R_earth): ", get_real(sim_param,"max_radius")/ExoplanetsSysSim.earth_radius)
    println(f, "# sigma_incl: ", get_real(sim_param,"sigma_incl"))
    println(f, "# sigma_incl_near_mmr: ", get_real(sim_param,"sigma_incl_near_mmr"))
    println(f, "# sigma_hk: ", get_real(sim_param,"sigma_hk"))
    println(f, "# num_mutual_hill_radii: ", get_real(sim_param,"num_mutual_hill_radii"))

    if string(get_function(sim_param,"generate_planet_mass_from_radius")) == "ExoplanetsSysSim.generate_planet_mass_from_radius_powerlaw"
        println(f, "# mr_power_index: ", get_real(sim_param,"mr_power_index"))
        println(f, "# mr_max_mass (M_earth): ", get_real(sim_param,"mr_max_mass")/ExoplanetsSysSim.earth_mass)
    elseif string(get_function(sim_param,"generate_planet_mass_from_radius")) == "generate_planet_mass_from_radius_Ning2018"
        println(f, "# mr_model: Ning2018")
    end

    println(f, "# sigma_log_radius_in_cluster: ", get_real(sim_param,"sigma_log_radius_in_cluster"))
    println(f, "# sigma_logperiod_per_pl_in_cluster: ", get_real(sim_param,"sigma_logperiod_per_pl_in_cluster"))
    println(f, "#")
end


