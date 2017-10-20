## ExoplanetsSysSim/examples/syssim_dist.jl
## (c) 2015 Eric B. Ford

# Calculaties the distance between simulation computed from an input "parameter file" and the summary statistics already stored in comp_file
# usage:  syssim_dist.jl param_file comp_file
# output: scalar_dist d_1 d_2 ...
#         where scalar_dist is the scalar distance between the simulated model's summary statistics and the precomputed summary statistics,
#         and d_1 d_2 ... represent a list of the individual distance metrics between the two sets of summary statistics

verbose = true

if length(ARGS)<2
   warn("usage:  syssim_dist.jl param_file comp_file")
   exit(255)
end

using Compat
#import Compat: UTF8String, ASCIIString

param_file = convert(String,ARGS[1])
if !isfile(param_file)
   warn(string("Can't read ",param_file))
   exit(255)
end

comp_file = convert(String,ARGS[2])
if !isfile(comp_file)
   warn(string("Can't read ",comp_file))
   exit(255)
end

if verbose println("# Loading ExoplanetsSysSim module."); end
using ExoplanetsSysSim

if verbose println("# Reading parameter file ", param_file,"."); end 
include(param_file)
if verbose println("# Calling setup_sim_param( ", ARGS[3:end], " )."); end
sim_param = setup_sim_param( convert(Array{String,1},ARGS[3:end]) )
if verbose println("# Active parameter values: ", make_vector_of_sim_param(sim_param) ); end

#add_param_fixed(sim_param,"rng_seed",1234)   # If you want to be able to reproduce simulations
if haskey(sim_param,"rng_seed")
   seed = get_int(sim_param,"rng_seed")
     if verbose println("# Seeding RNG with seed = ", seed); end
   srand(seed)
end

if verbose println("# Loading summary statistics from ", comp_file,"."); end
ss_comp = load_summary_stats(comp_file)

num_repeats = 1
for k in 1:num_repeats
if num_repeats > 1  print("\n\n"); end
if verbose println("# Generating physical catalog."); end
@time cat_phys = generate_kepler_physical_catalog(sim_param)
@time cat_phys = generate_kepler_physical_catalog(sim_param)

if verbose println("# Generating observational catalog."); end
@time cat_obs = observe_kepler_targets_single_obs(cat_phys,sim_param)
@time cat_obs = observe_kepler_targets_single_obs(cat_phys,sim_param)

ss = CatalogSummaryStatistics()
if haskey(sim_param,"use obs data only") &&  convert(Bool,get(sim_param,"use obs data only"))
  if verbose println("# Calculating summary statistics based on all data."); end
  @time ss = calc_summary_stats_obs_demo(cat_obs,sim_param)
else
  if verbose println("# Calculating summary statistics based on observational data only."); end
  @time ss = calc_summary_stats_sim_pass_one_demo(cat_obs,cat_phys,sim_param)
  @time ss = calc_summary_stats_sim_pass_two_demo(cat_obs,cat_phys,ss,sim_param)
end

if verbose println("# Calculating distance (pass 1)"); end
@time d1 = calc_distance_vector_demo(ss_comp, ss, 1, sim_param)
if verbose println("# Calculating distance (pass 2)"); end
@time d2 = calc_distance_vector_demo(ss_comp, ss, 2, sim_param)

d = [d1; d2]
@time dist = calc_scalar_distance(d)
print(dist, "   ")
print(STDOUT,join([d1;d2],' '))
print("\n")
end
