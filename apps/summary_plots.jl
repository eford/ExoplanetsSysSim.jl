# ExoplanetsSysSim/src/summary_plots.jl
## (c) 2015 Darin Ragozzine

try
if !(Pkg.installed("PyPlot") >= v"2.1.0")  # version chosen arbitrarily
  error("You need to upgrade the Optim package for this example.")
end
catch
  error("You need to install the Optim package for this example.")
end

verbose = true

export plot_all_summary_stats

if length(ARGS)<2
   warn("usage:  syssim_plots.jl param_file comp_file")
   println("ASSUMING!!: param_file = demo_param.in AND comp_file = demo_ss.out")

   ## This hack allows this command to be run from inside the REPL or at the 
   ## command line, though Julia understandably doesn't like it. 
   ARGS=["demo_param.in","demo_ss.out"]
   # exit(255)
end

param_file = convert(String,ARGS[1])
if !isreadable(param_file) || !isfile(param_file)
   warn(string("Can't read ",param_file))
   exit(255)
end

comp_file = convert(String,ARGS[2])
if !isreadable(comp_file) || !isfile(comp_file)
   warn(string("Can't read ",comp_file))
   exit(255)
end

if verbose println("# Loading ExoplanetsSysSim module."); end
using ExoplanetsSysSim
if verbose println("# Loading PyPlot Module; that's how we'll do plots."); end
using PyPlot
pygui(false)  # Tell PyPlot not to use the GUI, but to return objects
 
if verbose println("# Reading parameter file ", param_file,"."); end 
include(param_file)

if verbose println("# Calling setup_sim_param( ", ARGS[3:end], " )."); end
sim_param = setup_sim_param( convert(Array{String,1},ARGS[3:end]) )
if verbose println("# Active parameter values: ", make_vector_of_sim_param(sim_param) ); end

if verbose println("# Loading summary statistics from ", comp_file,"."); end
ss_comp = load_summary_stats(comp_file)

# Maybe? 
type PlotSummaryStatistics
  plotnames::Dict{String,Any}  	# For storing names of plots
  plotfilenames::Dict{String,Any}  # For storing filenames of plots
end
function PlotSummaryStatistics()
  PlotSummaryStatistics( Dict{String,Any}(), Dict{String,Any}() )
end

function plot_all_summary_stats(ss::CatalogSummaryStatistics, 
      sim_param::SimParam, verbose::Bool = false)
## TODO: Figure out the best filename to output to (user defined?) and implement.

 
# Print P vs. depth log-log scatter plot
if haskey(ss.stat,"P list") && haskey(ss.stat,"depth list")
	figure()
	## Some searching and outputting directly to a file without popping
	## up a figure window is not obvious, at least using PyPlot. Could use Gadfly.
	## TODO: Allow to output directly to a file without popping up a window
	loglog(ss.stat["P list"],ss.stat["depth list"],".")
    savefig("P_vs_depth.png")
else
	if verbose println("Cant plot P vs. depth because these are not found."); end
end #plot P vs. depth

# Plot multiplicity histogram in semilogy histogram-y scatter plot
if haskey(ss.stat,"num_sys_tranets")
	figure()
	semilogy(ss.stat["num_sys_tranets"],drawstyle="steps")
    savefig("num_sys_tranets.png")
else
	if verbose println("Cant plot multiplicity histogram because not found."); end
end # plot multiplicity histogram


end # plot_all_summary_stats


plot_all_summary_stats(ss_comp, sim_param, verbose)
