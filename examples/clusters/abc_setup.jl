## ExoplanetsSysSim/examples/setup_abc_christiansen.jl
## (c) 2016 Eric B. Ford & Danley C. Hsu

module EvalSysSimModel
  export setup, get_param_vector, get_ss_obs #, evaluate_model
  export gen_data, calc_summary_stats, calc_distance, is_valid
  using ExoplanetsSysSim
  include(joinpath(Pkg.dir(),"ExoplanetsSysSim","examples","multiplanet_systems", "model.jl"))
  include(joinpath(pwd(), "param_file.jl"))

  sim_param_closure = SimParam()
  summary_stat_ref_closure =  CatalogSummaryStatistics()
  dist_weights = ones(Float64,0)

    function is_valid(param_vector::Vector{Float64})
      global sim_param_closure
      update_sim_param_from_vector!(param_vector,sim_param_closure)
      if haskey(sim_param_closure,"mean_log_planet_radius")
        const mu_log_r::Float64 = get_real(sim_param_closure,"mean_log_planet_radius")
        const sigma_log_r::Float64 = get_real(sim_param_closure,"sigma_log_planet_radius")
        const min_radius::Float64 = get_real(sim_param_closure,"min_radius")
        const max_radius::Float64 = get_real(sim_param_closure,"max_radius")
        if sigma_log_r <= 0. || !(log(min_radius)<=mu_log_r<=log(max_radius)) 
         return false
        end
      end
      if haskey(sim_param_closure,"mean_log_planet_period")
        const mu_log_P::Float64 = get_real(sim_param_closure,"mean_log_planet_period")
        const sigma_log_P::Float64 = get_real(sim_param_closure,"sigma_log_planet_period")
        const min_period::Float64 = get_real(sim_param_closure,"min_period")
        const max_period::Float64 = get_real(sim_param_closure,"max_period")
        if sigma_log_P<=0. || !(log(min_period)<=mu_log_P<=log(max_period)) 
         return false
        end
      end
      if haskey(sim_param_closure,"frac_zero_planet")
         const frac_zero_planet::Float64 = get_real(sim_param_closure,"frac_zero_planet")
         if !(0.0<=frac_zero_planet<=1.0)
            return false
         end
      end
      if haskey(sim_param_closure,"frac_one_planet")
         const frac_one_planet::Float64 = get_real(sim_param_closure,"frac_one_planet")
         if !(0.0<=frac_one_planet<=1.0)
            return false
         end
      end

      if haskey(sim_param_closure,"fracs_num_planets")
         fracs::Array{Float64,1} = get(sim_param_closure,"fracs_num_planets",ones(1))
         if any(.!(0.0.<=fracs.<=1))
            return false
         end
      end

      if haskey(sim_param_closure,"eta_pl")
         eta_pl = get_real(sim_param_closure,"eta_pl")
         if any(!(0.0.<=eta_pl))
            return false
         end
      end

#      const rate_tab::Array{Float64,2} = get_any(sim_param_closure, "obs_par", Array{Float64,2})
#      const lambda = sum_kbn(rate_tab)
    # println("Lambda = ", lambda)
#      if sigma_log_r <= 0. || sigma_log_P<=0. || !(log(min_radius)<=mu_log_r<=log(min_radius)) || !(log(min_period)<=mu_log_P<=log(min_period)) || !(0.0<=frac_zero_planet<=1.0) || !(0.0<=frac_one_planet<=1.0) # || lambda > 10. || any(x -> x < 0., rate_tab)
#         return false
#      end
      return true
    end

#   Force the fracts_num_planets vector to be a non-negative simplex.  
#   If negative values happen often, then should change or adjust IS density to compensate
#   TODO FUTURE: I think Jessi has a better way to draw such vectors, but implementing that that might complicate the proposal step
    function normalize_categorical(param_vector::Vector{Float64}) 
      global sim_param_closure
      r::Range = get_range_for_sim_param("fracs_num_planets",sim_param_closure)
      param_vector[r] = abs(param_vector[r])   
      param_vector[r] = param_vector[r]/sum(param_vector[r])  
      update_sim_param_from_vector!(param_vector,sim_param_closure)

      if is_active(sim_param_closure,"sigma_incl")
        sig = get_real(sim_param_closure,"sigma_incl")
        if sig<0
           update_param(sim_param_closure,"sigma_incl",abs(sig))
        end
      end
      param_vector[1:end] = get_param_vector()
      return param_vector
    end

    function normalize_poisson_mixture(param_vector::Vector{Float64}) 
      global sim_param_closure
      update_sim_param_from_vector!(param_vector,sim_param_closure)
      frac_zero::Float64 = min(max(get_real(sim_param_closure,"frac_zero_planet"), 0.0),1.0)
      frac_one::Float64  = min(max(get_real(sim_param_closure,"frac_one_planet"), 0.0),1.0)
      frac_multi::Float64  = min(max(get_real(sim_param_closure,"frac_multi_planet"),0.0), 1.0)
      frac_sum = frac_zero+frac_one+frac_multi
      frac_zero /= frac_sum
      frac_one /= frac_sum
      frac_multi /= frac_sum
      update_param(sim_param_closure, "frac_zero_planet",frac_zero)
      update_param(sim_param_closure, "frac_one_planet",frac_one)
      update_param(sim_param_closure, "frac_multi_planet",frac_multi)

      if is_active(sim_param_closure,"sigma_incl")
        sig = get_real(sim_param_closure,"sigma_incl")
        if sig<0
           update_param(sim_param_closure,"sigma_incl",abs(sig))
        end
      end
      if is_active(sim_param_closure,"eta_pl")
        eta = get_real(sim_param_closure,"eta_pl")
        if eta<0
           update_param(sim_param_closure,"eta_pl",abs(eta))
        end
      end
      if is_active(sim_param_closure,"sigma_hk")
        sig = get_real(sim_param_closure,"sigma_hk")
        if sig<0
           update_param(sim_param_closure,"sigma_hk",abs(sig))
        end
      end
      param_vector[1:end] = get_param_vector()
      return param_vector
    end

#=
    function norm_log_simplex(param_vector::Vector{Float64})
      log_z_sum = logsumexp(abs(param_vector))
      param_vector = exp(abs(param_vector)-z_sum)
    end
=#

    function gen_data(param_vector::Vector{Float64})
      global sim_param_closure
      update_sim_param_from_vector!(param_vector,sim_param_closure)
      cat_phys = generate_kepler_physical_catalog(sim_param_closure)
      cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys, sim_param_closure)
      cat_obs = ExoplanetsSysSim.observe_kepler_targets_single_obs(cat_phys_cut, sim_param_closure)
      return cat_obs
    end

    function calc_summary_stats(cat::KeplerObsCatalog )
      global sim_param_closure
      sum_stat = calc_summary_stats_model(cat, sim_param_closure)
      return sum_stat
    end

    function calc_distance_no_weights(sum_stat_obs::CatalogSummaryStatistics,sum_stat_sim::CatalogSummaryStatistics, n::Integer = 0)
      global sim_param_closure
      dist = calc_distance_model(sum_stat_obs,sum_stat_sim, sim_param_closure)
      num_available = length(dist)
      num_to_use = n>0 ? min(n,num_available) : num_available
      dist_scalar = calc_scalar_distance(dist[1:num_to_use])
      println("param = ", get_param_vector()," dist = ",dist_scalar," : ", dist )
      return dist_scalar
    end
    function calc_distance(sum_stat_obs::CatalogSummaryStatistics,sum_stat_sim::CatalogSummaryStatistics, n::Integer = 0)
      global sim_param_closure, dist_weights
      dist = calc_distance_model(sum_stat_obs,sum_stat_sim, sim_param_closure)
      num_available = length(dist)
      num_to_use = n>0 ? min(n,num_available) : num_available
      @assert length(dist_weights)>=num_to_use
      idx_to_use = dist_weights[1:num_to_use] .> 0.0 
      dist_scalar = calc_scalar_distance(dist[idx_to_use].*dist_weights[idx_to_use])
      println("param = ", get_param_vector()," dist = ",dist_scalar," : ", dist, "; weighted = ", dist[idx_to_use].*dist_weights[idx_to_use])
      return dist_scalar
    end

  function calc_distance_scales(param::Vector{Float64}; n::Integer = 20)
     global sim_param_closure
     println("# param =",param)
     cat_obs = gen_data(param)
     sum_stat_obs = calc_summary_stats(cat_obs)
     dist_test = calc_distance_model(sum_stat_obs,sum_stat_obs,sim_param_closure)
     dist_arr = Array{Float64}(length(dist_test),n)
     for i in 1:n
        cat_sim = gen_data(param)
        sum_stat_sim = calc_summary_stats(cat_sim)
        dist_arr[:,i] = calc_distance_model(sum_stat_obs,sum_stat_sim,sim_param_closure)     
     end
     dist_medians = median(dist_arr,2)
     println("# Distance Means  =", mean(dist_arr,2))
     println("# Distance Medians=", median(dist_arr,2))
     println("# Distance Stddevs=", std(dist_arr,2))
     println("# 25th Percentile= ", map(i->quantile(dist_arr[i,:],0.25),1:size(dist_arr,1)))
     println("# 75th Percentile= ", map(i->quantile(dist_arr[i,:],0.75),1:size(dist_arr,1)))
     global dist_weights = zeros(Float64,length(dist_medians))
     for i in 1:length(dist_medians)
         if dist_medians[i]>0.0
            dist_weights[i] = 1.0/dist_medians[i]
         end       
     end
     println("# Distance Weights=", dist_weights)
  end

  function setup()
    global sim_param_closure = setup_sim_param_model()

    sim_param_closure = set_test_param(sim_param_closure)
    add_param_fixed(sim_param_closure,"stellar_catalog","q1_q17_dr25_stellar.jld") # "q1_q12_christiansen.jld")
    
    ### Use simulated planet candidate catalog data
    add_param_fixed(sim_param_closure,"num_kepler_targets",150969)  # For "observed" catalog
    #cat_obs = simulated_read_kepler_observations(sim_param_closure)
    #inv_det_prob(cat_obs, sim_param_closure)
    #println("------------------------")
    #global summary_stat_ref_closure = calc_summary_stats_model(cat_obs,sim_param_closure, trueobs_cat = true)
    ###
    
    ### Use real planet candidate catalog data
    add_param_fixed(sim_param_closure,"koi_catalog","q1_q16_koi_cand.csv")
    cat_obs = setup_actual_planet_candidate_catalog(setup_star_table_christiansen(sim_param_closure), sim_param_closure)
    global summary_stat_ref_closure = calc_summary_stats_model(cat_obs,sim_param_closure, trueobs_cat = true)
    ###

    num_targ_sim = floor(Int64,get_int(sim_param_closure,"num_kepler_targets")/10)  # TODO WARNING TEMPORARY: This is just for testing to speed things along...
    add_param_fixed(sim_param_closure,"num_targets_sim_pass_one",num_targ_sim)  # Set universal simulated catalog size

    calc_distance_scales(get_param_vector())
  end

  get_param_vector() = make_vector_of_sim_param(sim_param_closure)
  get_ss_obs() = summary_stat_ref_closure

  function set_simparam_ss(sim_param::ExoplanetsSysSim.SimParam, ss_true::ExoplanetsSysSim.CatalogSummaryStatistics)
    global sim_param_closure = sim_param
    global summary_stat_ref_closure = ss_true
  end


  function profile_model_eval(n::Integer)
    global sim_param_closure
    setup()
    global cat_phys_try_closure = generate_kepler_physical_catalog(sim_param_closure)
    cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys, sim_param_closure)
    global cat_obs_try_closure = observe_kepler_targets_single_obs(cat_phys_try_closure,sim_param_closure)
    global summary_stat_ref_closure = calc_summary_stats(cat_obs_try_closure)

    param_base = make_vector_of_sim_param(sim_param_closure)

    evaluate_model_demo_scalar_param_scalar_ret(param::Float64) = sum(evaluate_model([param;param_base[2:end]]))

    io = open("profile_syssim.tmp","w")
    Profile.print(io) # ,cols=800)
    close(io)
    #rm(io)

    Profile.init(10000000,0.001)
    eta_min = 0.5*exp(get_real(sim_param_closure,"log_eta_pl"))
    eta_max = 2.0*exp(get_real(sim_param_closure,"log_eta_pl"))
    #gc(); if VERSION >= VersionNumber(0,4,0)     gc_enable(false);   else     gc_disable();   end
    tic()
    for (i,eta) in enumerate(linspace(eta_min,eta_max,n) )
      local dist
      if i==2
         Profile.clear()
         Profile.clear_malloc_data()
      end
      @profile dist = evaluate_model_demo_scalar_param_scalar_ret(log(eta))
      println("eta= ",eta," d= ", dist)
    end
    toc()
    #gc_enable()
    pwd()
    io = open("profile_syssim.out","w")
    Profile.print(io,cols=800)
    close(io)
  end

end  # module EvalSysSimModel

include(joinpath(Pkg.dir("ABC"),"src/composite.jl"))

module SysSimABC
  export setup_abc, run_abc
  import ABC
  import Distributions
  using CompositeDistributions
  import ExoplanetsSysSim
  import EvalSysSimModel
  include(joinpath(Pkg.dir(),"ExoplanetsSysSim","examples","multiplanet_systems", "model.jl"))
  #include(joinpath(pwd(), "param_file.jl"))

  function setup_abc(num_dist::Integer = 0)
    EvalSysSimModel.setup()
    #println("# setup_abc: ",EvalSysSimModel.sim_param_closure)
    theta_true = EvalSysSimModel.get_param_vector()
    prior_dict = Dict("log_eta_pl" => Distributions.Uniform(0.1, log(10.)),"eta_pl" => Distributions.Uniform(0., 10.), "frac_zero_planet" => Distributions.Uniform(0.0, 1.0), "frac_one_planet" => Distributions.Uniform(0.0, 1.0), "frac_multi_planet" => Distributions.Uniform(0.0,1.0),
       "mean_log_planet_radius" => Distributions.Uniform(log(0.5*ExoplanetsSysSim.earth_radius), log(4.*ExoplanetsSysSim.earth_radius)), "mean_log_planet_period" => Distributions.Uniform(log(1.), log(100.)), 
       "sigma_log_planet_radius" => Distributions.Uniform(0., log(4.)), "sigma_log_planet_period" => Distributions.Uniform(0., log(4.)), 
       "power_law_P" => Distributions.Uniform(-1.0, 1.0), "power_law_r" => Distributions.Uniform(-4.0, -1.0),
       "sigma_hk" => Distributions.Uniform(0.0,1.0),  "sigma_hk_one" => Distributions.Uniform(0.0,1.0),  "sigma_hk_multi" => Distributions.Uniform(0.0,1.0),
       "sigma_incl" => Distributions.Uniform(0.0,30.0)   # degrees; 0 = coplanar w/ generate_kepler_target_simple; ignored by generate_planetary_system_uncorrelated_incl
)
    active_param_keys_sorted = make_vector_of_active_param_keys(EvalSysSimModel.sim_param_closure)
    array_of_priors = Array{Distributions.ContinuousDistribution}(length(active_param_keys_sorted) )
    for i in 1:length(array_of_priors) 
        array_of_priors[i]  = prior_dict[active_param_keys_sorted[i]];
    end
    param_prior = CompositeDist( array_of_priors ) 
    # param_prior = CompositeDist( Distributions.ContinuousDistribution[Distributions.Uniform(0., 2.) for x in 1:length(theta_true)] )
    in_parallel = nworkers() > 1 ? true : false

    calc_distance_ltd(sum_stat_obs::ExoplanetsSysSim.CatalogSummaryStatistics,sum_stat_sim::ExoplanetsSysSim.CatalogSummaryStatistics) = EvalSysSimModel.calc_distance(sum_stat_obs,sum_stat_sim,num_dist)

    global abc_plan = ABC.abc_pmc_plan_type(EvalSysSimModel.gen_data,EvalSysSimModel.calc_summary_stats,calc_distance_ltd, param_prior, make_proposal_dist=ABC.make_proposal_dist_gaussian_diag_covar, is_valid=EvalSysSimModel.is_valid,
         #normalize=EvalSysSimModel.normalize_poisson_mixture,
         normalize=EvalSysSimModel.normalize_categorical,
                                     num_part=20, num_max_attempt=200, num_max_times=60, epsilon_init=9.9e99, target_epsilon=1.0e-100, in_parallel=in_parallel, adaptive_quantiles = true, epsilon_reduction_factor=0.8, tau_factor=2.0);
  end

  function run_abc_largegen(pop::ABC.abc_population_type, ss_true::ExoplanetsSysSim.CatalogSummaryStatistics, epshist_targ::Float64, npart::Integer = 1000, num_dist::Integer = 0)
    sim_param_closure = setup_sim_param_model()
    sim_param_closure = EvalSysSimModel.set_test_param(sim_param_closure)
    EvalSysSimModel.set_simparam_ss(sim_param_closure, ss_true)	

    theta_true = EvalSysSimModel.get_param_vector()
    prior_dict = Dict("log_eta_pl" => Distributions.Uniform(0.1, log(10.)), "eta_pl" => Distributions.Uniform(0., 10.),"frac_zero_planet" => Distributions.Uniform(0.0, 1.0), "frac_one_planet" => Distributions.Uniform(0.0, 1.0), "frac_multi_planet" => Distributions.Uniform(0.0,1.0),
       "fracs_num_planets" => Distributions.Dirichlet(10,1.0),
       "mean_log_planet_radius" => Distributions.Uniform(log(0.5*ExoplanetsSysSim.earth_radius), log(4.*ExoplanetsSysSim.earth_radius)), "mean_log_planet_period" => Distributions.Uniform(log(1.), log(100.)), 
       "sigma_log_planet_radius" => Distributions.Uniform(0., log(4.)), "sigma_log_planet_period" => Distributions.Uniform(0., log(4.)), 
       "power_law_P" => Distributions.Uniform(-1.0, 1.0), "power_law_r" => Distributions.Uniform(-4.0, -1.0),
       "sigma_hk" => Distributions.Uniform(0.0,1.0),  "sigma_hk_one" => Distributions.Uniform(0.0,1.0),  "sigma_hk_multi" => Distributions.Uniform(0.0,1.0),
)
    active_param_keys_sorted = make_vector_of_active_param_keys(EvalSysSimModel.sim_param_closure)
    array_of_priors = Array{Distributions.ContinuousDistribution}(length(active_param_keys_sorted) )
    for i in 1:length(array_of_priors) 
        array_of_priors[i]  = prior_dict[active_param_keys_sorted[i]];
    end
    param_prior = CompositeDist( array_of_priors ) 
    # param_prior = CompositeDist( Distributions.ContinuousDistribution[Distributions.Uniform(0., 2.) for x in 1:length(theta_true)] )
    in_parallel = nworkers() > 1 ? true : false

    calc_distance_ltd(sum_stat_obs::ExoplanetsSysSim.CatalogSummaryStatistics,sum_stat_sim::ExoplanetsSysSim.CatalogSummaryStatistics) = EvalSysSimModel.calc_distance(sum_stat_obs,sum_stat_sim,num_dist)

    global abc_plan = ABC.abc_pmc_plan_type(EvalSysSimModel.gen_data,EvalSysSimModel.calc_summary_stats,calc_distance_ltd, param_prior, is_valid=EvalSysSimModel.is_valid,
                                     num_part=npart, num_max_attempt=200, num_max_times=1, epsilon_init=9.9e99, target_epsilon=1.0e-100, in_parallel=in_parallel);

    println("# run_abc_largegen: ",EvalSysSimModel.sim_param_closure)
    sampler_plot = abc_plan.make_proposal_dist(pop, abc_plan.tau_factor)
    theta_plot = Array{Float64}(size(pop.theta, 1), npart)
    for i in 1:npart
      theta_plot[:,i], dist_plot, attempts_plot = ABC.generate_theta(abc_plan, sampler_plot, ss_true, epshist_targ)
    end
    return theta_plot
  end

  function run_abc(abc_plan::ABC.abc_pmc_plan_type)
    println("# run_abc: ",EvalSysSimModel.sim_param_closure)
    ss_true = EvalSysSimModel.get_ss_obs()
    #println("True catalog SS: ", ss_true)
    pop_out = ABC.run_abc(abc_plan,ss_true;verbose=true)
  end

  function run_abc(abc_plan::ABC.abc_pmc_plan_type, pop::ABC.abc_population_type)
    dist_threshold = maximum(pop.dist)
    EvalSysSimModel.add_param_fixed(EvalSysSimModel.sim_param_closure,"minimum ABC dist skip pass 2",dist_threshold)
    println("# run_abc: ",EvalSysSimModel.sim_param_closure)
    ss_true = EvalSysSimModel.get_ss_obs()
    pop_out = ABC.run_abc(abc_plan,ss_true,pop;verbose=true)
  end

end # module SysSimABC
