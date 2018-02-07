## ExoplanetsSysSim/examples/profile.jl
## (c) 2015 Eric B. Ford

#workspace()
if !isdefined(:ExoplanetsSysSim)
  using ExoplanetsSysSim
end

module EvalModel
  export profile_model_eval
  using ExoplanetsSysSim

  sim_param_closure = SimParam()
  cat_phys_try_closure = KeplerPhysicalCatalog([])
  cat_obs_try_closure = KeplerObsCatalog([])
  summary_stat_try_closure =  CatalogSummaryStatistics()
  summary_stat_ref_closure =  CatalogSummaryStatistics()

  include(joinpath(Pkg.dir("ExoplanetsSysSim"), "src", "eval_model_macro.jl"))
  @make_evaluate_model(sim_param_closure,cat_phys_try_closure,cat_obs_try_closure,summary_stat_try_closure,summary_stat_ref_closure)

  function profile_model_eval(n::Integer)
    global sim_param_closure = setup_sim_param_demo()
    add_param_fixed(sim_param_closure,"eta_pl",0.3)
    global cat_phys_try_closure = generate_kepler_physical_catalog(sim_param_closure)
    global cat_obs_try_closure = observe_kepler_targets_single_obs(cat_phys_try_closure,sim_param_closure)
    global summary_stat_ref_closure = calc_summary_stats_obs_demo(cat_obs_try_closure,sim_param_closure)

    param_base = make_vector_of_sim_param(sim_param_closure)

    evaluate_model_demo_scalar_param_scalar_ret(param::Float64) = sum(evaluate_model([param;param_base[2:end]]))

    io = open(joinpath(Pkg.dir("ExoplanetsSysSim"), "examples", "profile_syssim.tmp"),"w")
    Profile.print(io) # ,cols=800)
    close(io)
    #rm(io)

    Profile.init(10000000,0.001)
    eta_min = 0.5*get_real(sim_param_closure,"eta_pl")
    eta_max = 2.0*get_real(sim_param_closure,"eta_pl")
    #gc(); if VERSION >= VersionNumber(0,4,0)     gc_enable(false);   else     gc_disable();   end
    tic()
    println("# eta_min = ", eta_min)
    println("# eta_max = ", eta_max)
    println("# enum = ", enumerate(linspace(eta_min,eta_max,n) ))
    for (i,eta) in enumerate(linspace(eta_min,eta_max,n) )
      local dist
      if i==2
         Profile.clear()
         Profile.clear_malloc_data()
      end
      @profile dist = evaluate_model_demo_scalar_param_scalar_ret(eta)
      println("eta= ",eta," d= ", dist)
    end
    toc()
    #gc_enable()
    pwd()
    io = open(joinpath(Pkg.dir("ExoplanetsSysSim"), "examples", "profile_syssim.out"),"w")
    Profile.print(io) # ,cols=800)
    close(io)
  end

end  # module EvalModel

using EvalModel

profile_model_eval(2)


