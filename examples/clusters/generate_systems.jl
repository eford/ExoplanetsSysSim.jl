include("clusters.jl")
sim_param = setup_sim_param_model()



tic()

##### To generate the underlying systems:

cat_phys = generate_kepler_physical_catalog(sim_param)

##### For saving the underlying/true planets/systems:

f = open("periods_all.out", "w")
write_model_params(f, sim_param)
for (i,targ) in enumerate(cat_phys.target)
    if length(targ.sys[1].orbit) > 0
        periods_sys = Array{Float64}(length(targ.sys[1].orbit))
        for (j,planet) in enumerate(targ.sys[1].orbit)
            periods_sys[j] = planet.P #days
        end
        println(f, periods_sys)
    end
end
close(f)

f = open("eccentricities_all.out", "w")
write_model_params(f, sim_param)
for (i,targ) in enumerate(cat_phys.target)
    if length(targ.sys[1].orbit) > 0
        ecc_sys = Array{Float64}(length(targ.sys[1].orbit))
        for (j,planet) in enumerate(targ.sys[1].orbit)
            ecc_sys[j] = planet.ecc
        end
        println(f, ecc_sys)
    end
end
close(f)

f = open("radii_all.out", "w")
write_model_params(f, sim_param)
for (i,targ) in enumerate(cat_phys.target)
    if length(targ.sys[1].planet) > 0
        radii_sys = Array{Float64}(length(targ.sys[1].planet))
        for (j,planet) in enumerate(targ.sys[1].planet)
            radii_sys[j] = planet.radius #solar radii
        end
        println(f, radii_sys)
    end
end
close(f)

f = open("masses_all.out", "w")
write_model_params(f, sim_param)
for (i,targ) in enumerate(cat_phys.target)
    if length(targ.sys[1].planet) > 0
        masses_sys = Array{Float64}(length(targ.sys[1].planet))
        for (j,planet) in enumerate(targ.sys[1].planet)
            masses_sys[j] = planet.mass #solar masses
        end
        println(f, masses_sys)
    end
end
close(f)



##### To generate the observed systems:

cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param)
cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param)
summary_stat = calc_summary_stats_model(cat_obs,sim_param)

toc()

##### For saving the observed planets/systems:

f = open("periods.out", "w")
write_model_params(f, sim_param)
for num_pl_in_sys in 1:length(summary_stat.cache["idx_n_tranets"])
  num_targets = length(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
  period_array = Array{Float64}(num_pl_in_sys,num_targets)
  for (i,j) in enumerate(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
      period_array[:,i] = map(ExoplanetsSysSim.period, cat_obs.target[j].obs)[1:num_pl_in_sys]
  end
  println(f,"# Periods of systems with ", num_pl_in_sys, " detected planets.")
  if length(period_array) > 0
    println(f,period_array')
  end
end
close(f)

f = open("depths.out", "w")
write_model_params(f, sim_param)
for num_pl_in_sys in 1:length(summary_stat.cache["idx_n_tranets"])
  num_targets = length(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
  depths_array = Array{Float64}(num_pl_in_sys,num_targets)
  for (i,j) in enumerate(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
      depths_array[:,i] = map(ExoplanetsSysSim.depth, cat_obs.target[j].obs)[1:num_pl_in_sys]
  end
  println(f,"# Transit depths of systems with ", num_pl_in_sys, " detected planets.")
  if length(depths_array) > 0
    println(f,depths_array')
  end
end
close(f)

f = open("durations.out", "w")
write_model_params(f, sim_param)
for num_pl_in_sys in 1:length(summary_stat.cache["idx_n_tranets"])
  num_targets = length(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
  duration_array = Array{Float64}(num_pl_in_sys,num_targets)
  for (i,j) in enumerate(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
      duration_array[:,i] = map(ExoplanetsSysSim.duration, cat_obs.target[j].obs)[1:num_pl_in_sys]
  end
  println(f,"# Transit durations of systems with ", num_pl_in_sys, " detected planets.")
  if length(duration_array) > 0
    println(f,duration_array')
  end
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
  if length(ecc_array) > 0
    println(f,ecc_array')
  end
end
close(f)
=#

