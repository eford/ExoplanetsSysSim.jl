## ExoplanetsSysSim/examples/setup_abc.jl
## (c) 2015 Eric B. Ford

# Add "minimum ABC dist skip pass 2"

module EvalSysSimModel
  export setup, get_param_vector, get_ss_obs #, evaluate_model
  export gen_data, calc_summary_stats, calc_distance, is_valid
  using ExoplanetsSysSim

  sim_param_closure = SimParam()
  summary_stat_ref_closure =  CatalogSummaryStatistics()

  # If we just wanted one call to evaluate_model, we could just use the macro
  #include("eval_model_macro.jl")
  #@make_evaluate_model(sim_param_closure,cat_phys_try_closure,cat_obs_try_closure,summary_stat_try_closure,summary_stat_ref_closure)
  # But we want a little finer control, so we'll do it manually....

    function is_valid(param_vector::Vector{Float64})
      global sim_param_closure
      update_sim_param_from_vector!(param_vector,sim_param_closure)
      const eta_pl = exp(get_real(sim_param_closure,"log_eta_pl"))
      const sigma_log_r = get_real(sim_param_closure,"sigma_log_planet_radius")
      const sigma_log_P = get_real(sim_param_closure,"sigma_log_planet_period")
      if   eta_pl <=0. || sigma_log_r <= 0. || sigma_log_P<=0.
         return false
      end
      return true
    end

    function gen_data(param_vector::Vector{Float64})
      global sim_param_closure
      update_sim_param_from_vector!(param_vector,sim_param_closure)
      cat_phys = generate_kepler_physical_catalog(sim_param_closure)
      cat_obs = observe_kepler_targets_single_obs(cat_phys,sim_param_closure)
      (cat_obs,cat_phys)
    end

    # TODO OPT: Eventually, could adapt ABC.jl to use distance from first pass to decide if should compute additional summary statistics
    function calc_summary_stats(cat::Tuple{KeplerObsCatalog,KeplerPhysicalCatalog} )
      global sim_param_closure
      # cat_obs = cat[1]
      # cat_phys = cat[2]
      sum_stat = calc_summary_stats_sim_pass_one_demo(cat[1],cat[2],sim_param_closure)
      sum_stat = calc_summary_stats_sim_pass_two_demo(cat[1],cat[2],sum_stat,sim_param_closure)
    end

    function calc_distance(sum_stat_obs::CatalogSummaryStatistics,sum_stat_sim::CatalogSummaryStatistics, n::Integer = 0)
      global sim_param_closure
      dist1 = calc_distance_vector_demo(sum_stat_obs,sum_stat_sim, 1, sim_param_closure)
      if haskey(sim_param_closure,"minimum ABC dist skip pass 2")
         if calc_scalar_distance(dist1) > get_real(sim_param_closure,"minimum ABC dist skip pass 2")
            return dist1
         end
      end
      dist2 = calc_distance_vector_demo(sum_stat_obs,sum_stat_sim, 2, sim_param_closure)
      num_available = length(dist1)+length(dist2)
      num_to_use = n>0 ? min(n,num_available) : num_available
      return calc_scalar_distance([dist1; dist2][1:num_to_use])
    end

  function setup()
    global sim_param_closure = setup_sim_param_demo()
    # To keep the abc demo runtime reasonable, use a modest number of target
    add_param_fixed(sim_param_closure,"num_targets_sim_pass_one",1900)  # For faster simulated catalogs
    add_param_fixed(sim_param_closure,"num_kepler_targets",190000)  # For "observed" catalog
    # Since we haven't parallelized the reading of stellar catalog files, the abc_demo uses made up targets, allowing us to run on multiple processors
    add_param_fixed(sim_param_closure,"generate_kepler_target",ExoplanetsSysSim.generate_kepler_target_simple)

    # cat_phys = generate_kepler_physical_catalog(sim_param_closure)
    # cat_obs = observe_kepler_targets_single_obs(cat_phys,sim_param_closure) 
    cat_obs = simulated_read_kepler_observations(sim_param_closure)
    global summary_stat_ref_closure = calc_summary_stats_obs_demo(cat_obs,sim_param_closure)
  end

  get_param_vector() = make_vector_of_sim_param(sim_param_closure)
  get_ss_obs() = summary_stat_ref_closure

end  # module EvalSysSimModel

module SysSimABC
  export setup_abc, run_abc
  import ABC
  import Distributions
  import ExoplanetsSysSim
  import EvalSysSimModel

  function setup_abc(num_dist::Integer = 0)
    global sim_param_closure
    EvalSysSimModel.setup()
    println("# setup_abc: ",EvalSysSimModel.sim_param_closure)
    theta_true = EvalSysSimModel.get_param_vector()
    param_prior = Distributions.MvNormal(theta_true,ones(length(theta_true)))
    in_parallel = nworkers() > 1 ? true : false

    calc_distance_ltd(sum_stat_obs::ExoplanetsSysSim.CatalogSummaryStatistics,sum_stat_sim::ExoplanetsSysSim.CatalogSummaryStatistics) = EvalSysSimModel.calc_distance(sum_stat_obs,sum_stat_sim,num_dist)


    abc_plan = ABC.abc_pmc_plan_type(EvalSysSimModel.gen_data,EvalSysSimModel.calc_summary_stats,calc_distance_ltd, param_prior, is_valid=EvalSysSimModel.is_valid,
                                     num_part=40, num_max_attempt=100, num_max_times=5, epsilon_init=9.9e99, target_epsilon=0.001, in_parallel=in_parallel);
  end

  function run_abc(abc_plan::ABC.abc_pmc_plan_type)
    #global sim_param_closure
    println("# run_abc: ",EvalSysSimModel.sim_param_closure)
    ss_true = EvalSysSimModel.get_ss_obs()
    pop_out = ABC.run_abc(abc_plan,ss_true;verbose=true)
  end

  function run_abc(abc_plan::ABC.abc_pmc_plan_type, pop::ABC.abc_population_type)
    #global sim_param_closure
    dist_threshold = maximum(pop.dist)
    ExoplanetsSysSim.add_param_fixed(EvalSysSimModel.sim_param_closure,"minimum ABC dist skip pass 2",dist_threshold)
    println("# run_abc: ",EvalSysSimModel.sim_param_closure)
    ss_true = EvalSysSimModel.get_ss_obs()
    pop_out = ABC.run_abc(abc_plan,ss_true,pop;verbose=true)
  end


end # module SysSimABC


