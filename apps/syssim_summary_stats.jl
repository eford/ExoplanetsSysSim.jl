## ExoplanetsSysSim/examples/syssim_summary_stats.jl
## (c) 2015 Eric B. Ford

# Calculaties summary statistics based on an input "parameter file" and outputs results to stdout or a file
# usage:  syssim_summary_stats.jl param_file [output_file]
# output: a binary (HDF5/JLD) file containing the summary statistics computed based on the parameter file.

verbose = true

if length(ARGS)<2
   warn("usage:  syssim_summary_stats.jl param_file output_file ")
   exit(255)
end

using Compat
import Compat: UTF8String, ASCIIString

param_file = convert(Compat.ASCIIString,ARGS[1])
if !isfile(param_file)
   warn(string("Can't read ",param_file))
   exit(255)
end

if verbose println("# Loading ExoplanetsSysSim module."); end
using ExoplanetsSysSim

if verbose println("# Reading parameter file ", param_file,"."); end 
include(param_file)
if verbose println("# Calling setup_sim_param( ", ARGS[3:end], " )."); end
sim_param = setup_sim_param( convert(Array{Compat.ASCIIString,1},ARGS[3:end]) )
if verbose println("# Active parameter values: ", make_vector_of_sim_param(sim_param) ); end

if haskey(sim_param,"rng_seed")
   seed = get_int(sim_param,"rng_seed")
     if verbose println("# Seeding RNG with seed = ", seed); end
   srand(seed)
end

if verbose println("# Generating physical catalog."); end
cat_phys = generate_kepler_physical_catalog(sim_param)

if verbose println("# Generating observational catalog."); end
cat_obs = observe_kepler_targets_single_obs(cat_phys,sim_param)

ss = CatalogSummaryStatistics()
if haskey(sim_param,"use obs data only") &&  convert(Bool,get(sim_param,"use obs data only"))
  if verbose println("# Calculating summary statistics based on all data."); end
  ss = calc_summary_stats_obs_demo(cat_obs,sim_param)
else
  if verbose println("# Calculating summary statistics based on observational data only."); end
  ss = calc_summary_stats_sim_pass_one_demo(cat_obs,cat_phys,sim_param)
  ss = calc_summary_stats_sim_pass_two_demo(cat_obs,cat_phys,ss,sim_param)
end

if length(ARGS)>=2 
  output_file = convert(Compat.ASCIIString,ARGS[2])
  if verbose println("# Writing summary statistics to ", output_file,"."); end
  try 
     save_sim_results(output_file,sim_param,summary_stats=ss)
  catch
     warn(string("Error writing to ",output_file))
  end
else  # Intetionally don't let you get here, since currently the "summary statistics" are so long, it's obnoxious to print.
      # in principle, could specialize how to print a summary statistics object in a useful way.  (I think the relevant function might be "show".)
  if verbose println("# Summary statistics:"); end
  println(ss)
end

