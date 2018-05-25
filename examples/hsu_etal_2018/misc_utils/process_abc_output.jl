##ExoplanetsSysSim/examples/hsu_etal_2018/misc_utils/process_abc_output.jl
##(c) 2017 Keir A. Ashby
# takes the output from abc_run_christiansen.jl and makes relevant plots 
# 1. make_plots sends the output into whatever plotting functions we choose to call
# 2. generation_plot makes a plot for the mean, median and overall distribution of distances or values of theta over each generation of abc. saves each plot

using PyPlot
using ExoplanetsSysSim
include(joinpath(Pkg.dir(), "ExoplanetsSysSim", "src", "constants.jl"))

#sends the output into two plotting functions
function make_plots(output::ABC.abc_population_type, ss_true::ExoplanetsSysSim.CatalogSummaryStatistics, run::Int64 = 1) #run indicates how many times abc_run_christiansen.jl has generated an output in the event that there are multiple outputs. 
  num_thetas = length(output.accept_log.theta[1]) 
#  gen_start = output.accept_log.generation_starts_at
  run = generation_plot(output.accept_log, output.accept_log.dist, run) 	#plots changing distances over each generation
  for i in 1:num_thetas
    run = generation_plot(output.accept_log, output.accept_log.theta, run, i) 	#plots changing theta values over each generation
  end 
end

#takes a log from the output and returns the mean, median, and overall distribution of either distance or a theta parameter over every generation of ABC
function generation_plot(log::ABC.abc_log_type, log_part::Array, run::Int64=1, m::Int64=1) #user can decide which part of theta she wants by specifying m in the input
  if typeof(log_part[1]) == Array{Float64,1} 					#first check whether the user has input an array of distances or of theta 
    my_theta = zeros(length(log_part))
    for i in 1:length(my_theta)
      my_theta[i] = log_part[i][m]
    end
    log_part = my_theta
    log_name = "Theta_$m";
    run +=10*m; 								#we keep different figures separated by factors of 10 so as not to overlap
  else
    log_name = "Distance";
  end

  avg_list = zeros(length(log.generation_starts_at)) 				#array holding the average of log_part over each generation of abc
  med_list = zeros(length(log.generation_starts_at)) 				#array holding the median of log_part over each generation of abc
  all_lengths = zeros(length(log.generation_starts_at)) 			#array holding the number of log_part values in each generation of abc
  for i in 1:length(log.generation_starts_at)   
    n = 1
    if i == length(log.generation_starts_at) 				
    gen_list = zeros(length(log_part)+1-log.generation_starts_at[i]) 		#array holding each value of log_part in a given generation
	for j in log.generation_starts_at[i]:length(log_part)
	  gen_list[n] = log_part[j]
	  n +=1
	end
    else    
    gen_list = zeros(log.generation_starts_at[i+1]-log.generation_starts_at[i]) #array holding each value of log_part in a given generation
	for j in log.generation_starts_at[i]:log.generation_starts_at[i+1]-1
	  gen_list[n] = log_part[j]
 	  n +=1
	end
    end
    avg_list[i] = mean(gen_list) 						#take the average of log_part over this generation
    med_list[i] = median(gen_list) 						#take the median of log_part over this generation
    all_lengths[i] = length(gen_list) 						#take the number of values of log_part for this generation
  end
  all_gen = ones(1) 								#array holding the generation number corresponding to every value of log_part
  for k in 1:length(log.generation_starts_at)
    this_gen = ones(convert(Int64, all_lengths[k])).*k
    append!(all_gen, this_gen)
  end
  shift!(all_gen)

  figure(run)
  semilogy(1:length(log.generation_starts_at), avg_list, "r") 			#plot of average log_part over each generation of ABC
  title(string(log_name, " vs. Generation"), fontsize=18)
  xlabel("Generation", fontsize=14)
  ylabel(log_name, fontsize=20)
  tick_params(labelsize=13, width=2)

  semilogy(1:length(log.generation_starts_at), med_list, "g") 			#plot of median log_part over each generation of ABC

  xmax= convert(Float64, length(log.generation_starts_at))+0.2 			#just so we can see everything
  semilogy(all_gen, log_part, ".")
  axis([0.8, xmax, 0, maximum(log_part)+0.01])

  savefig(string("plot_", log_name, ".png")) 
  return(run)
end
