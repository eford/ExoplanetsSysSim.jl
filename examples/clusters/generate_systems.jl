include("clusters.jl")
sim_param = setup_sim_param_model()
cat_phys = generate_kepler_physical_catalog(sim_param)
cat_obs = observe_kepler_targets_single_obs(cat_phys,sim_param)
summary_stat = calc_summary_stats_model(cat_obs,sim_param)

f = open("periods.out", "w")
for num_pl_in_sys in 1:length(summary_stat.cache["idx_n_tranets"])
  num_targets = length(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
  period_array = Array{Float64}(num_pl_in_sys,num_targets)
  for (i,j) in enumerate(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
      period_array[:,i] = map(ExoplanetsSysSim.period, cat_obs.target[j].obs)[1:num_pl_in_sys]
  end
  println(f,"# Periods of systems with ", num_pl_in_sys, " detected planets.")
  println(f,period_array')
end
close(f)

f = open("depths.out", "w")
for num_pl_in_sys in 1:length(summary_stat.cache["idx_n_tranets"])
  num_targets = length(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
  depths_array = Array{Float64}(num_pl_in_sys,num_targets)
  for (i,j) in enumerate(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
      depths_array[:,i] = map(ExoplanetsSysSim.depth, cat_obs.target[j].obs)[1:num_pl_in_sys]
  end
  println(f,"# Transit depths of systems with ", num_pl_in_sys, " detected planets.")
  println(f,depths_array')
end
close(f)

f = open("durations.out", "w")
for num_pl_in_sys in 1:length(summary_stat.cache["idx_n_tranets"])
  num_targets = length(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
  duration_array = Array{Float64}(num_pl_in_sys,num_targets)
  for (i,j) in enumerate(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
      duration_array[:,i] = map(ExoplanetsSysSim.duration, cat_obs.target[j].obs)[1:num_pl_in_sys]
  end
  println(f,"# Transit durations of systems with ", num_pl_in_sys, " detected planets.")
  println(f,duration_array')
end
close(f)

#= 
f = open("eccentricities.out", "w")
for num_pl_in_sys in 1:length(summary_stat.cache["idx_n_tranets"])
  num_targets = length(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
  ecc_array = Array{Float64}(num_pl_in_sys,num_targets)
  for (i,j) in enumerate(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
      ecc_array[:,i] = map(ExoplanetsSysSim.eccentricity, cat_obs.target[j].obs)[1:num_pl_in_sys]
  end
  println(f,"# Eccentricities of systems with ", num_pl_in_sys, " detected planets.")
  println(f,ecc_array')
end
close(f)
=#

