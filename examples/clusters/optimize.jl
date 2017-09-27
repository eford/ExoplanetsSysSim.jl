include("clusters.jl")
sim_param = setup_sim_param_model()
add_param_fixed(sim_param,"num_targets_sim_pass_one",150969)   # For "observed" data, use a realistic number of targets (after any cuts you want to perform)

# For now we're creating a simualted catalog that we fit to.  
# One of Danley/Darin/Keir probably have code to read in a Kepler catalog, so can compute summary statistics from it.
# Or for now, you could just manually set arrays with the properties to be compared in calc_distance

cat_phys = generate_kepler_physical_catalog(sim_param)
cat_obs = observe_kepler_targets_single_obs(cat_phys,sim_param)
summary_stat_ref = calc_summary_stats_model(cat_obs,sim_param)


#add_param_fixed(sim_param,"num_targets_sim_pass_one",150969*2)   # Could change the number of targets used for simulated catalogs if that helps by reducing Monte Carlo variance


function calc_distance(ss1::ExoplanetsSysSim.CatalogSummaryStatistics, ss2::ExoplanetsSysSim.CatalogSummaryStatistics)
  d = Array{Float64}(6)
  d[1] = abs( ss1.stat["num_tranets"]/ss1.stat["num targets"] - ss2.stat["num_tranets"]/ss2.stat["num targets"] ) 
  d[2] = ksstats(ss1.stat["P list"],ss2.stat["P list"])[5]
  d[3] = ksstats(ss1.stat["depth list"],ss2.stat["depth list"])[5] 
  d[4] = ksstats(ss1.stat["num n-tranet systems"],ss2.stat["num n-tranet systems"])[5]
  d[5] = ksstats(ss1.stat["period_ratio_list"],ss2.stat["period_ratio_list"])[5]
  d[6] = ksstats(ss1.stat["duration_ratio_list"],ss2.stat["duration_ratio_list"])[5]
  return sum(d)
end

function target_function(active_param::Vector)
  global sim_param, summary_stat_ref
  sim_param_here = deepcopy(sim_param)
  ExoplanetsSysSim.update_sim_param_from_vector!(active_param,sim_param_here)
  cat_phys = generate_kepler_physical_catalog(sim_param_here)
  cat_obs = observe_kepler_targets_single_obs(cat_phys,sim_param_here)
  summary_stat = calc_summary_stats_model(cat_obs,sim_param_here)

  dist = calc_distance(summary_stat,summary_stat_ref)
end


tic()
active_param_true = make_vector_of_sim_param(sim_param)
println("# True values: ", active_param_true)
num_eval = 20
results = map(x->target_function(active_param_true), 1:num_eval)
mean_dist = mean(results)
rms_dist = std(results)
println("# Distance using true values: ", mean_dist, " +/- ",rms_dist)
toc()

active_param = 2*active_param_true
#=  
using Optim   # I didn't have such good results with my initial attempts with BFGS, so I moved on
opt_result = optimize(target_function, active_param, method=BFGS(), f_tol=rms_dist, allow_f_increases=true, show_trace=true)
=#

# Pkg.add("BlackBoxOptim")       # only need to do these once
# Pkg.checkout("BlackBoxOptim")  # needed to get the lastest version
using BlackBoxOptim              # see https://github.com/robertfeldt/BlackBoxOptim.jl for documentation

# Eventually, we should make these physically motivated limits when we don't know the true values
sr = [(active_param_true[i]/2,active_param_true[i]*2) for i in 1:length(active_param)]  

# May want to experiment with different algorithms, number of evaluations, population size, etc.
tic()
opt_result = bboptimize(target_function; SearchRange = sr, NumDimensions = length(active_param), Method = :adaptive_de_rand_1_bin_radiuslimited, PopulationSize = length(active_param)*4, MaxFuncEvals = 100, TraceMode = :verbose)  

#opt_result = bboptimize(target_function; SearchRange = sr, NumDimensions = length(active_param), Method = :adaptive_de_rand_1_bin_radiuslimited, PopulationSize = length(active_param)*4, MaxFuncEvals = 100, TargetFitness = mean_dist, FitnessTolerance = rms_dist/10, TraceMode = :verbose)  
toc()

