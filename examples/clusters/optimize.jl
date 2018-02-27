import DataArrays.skipmissing

include("clusters.jl")

sim_param = setup_sim_param_model()
add_param_fixed(sim_param,"num_targets_sim_pass_one",15006) #9)   # For "observed" data, use a realistic number of targets (after any cuts you want to perform)

##### To generate a simulated catalog to fit to:

cat_phys = generate_kepler_physical_catalog(sim_param)
cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param)
cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param)
summary_stat_ref = calc_summary_stats_model(cat_obs,sim_param)


##### To load the DR25 catalog and compute arrays of multiplicities, periods, period ratios, transit durations, transit depths, period-normalized transit duration ratios (xi), and transit depth ratios:

Q1Q17_DR25 = CSV.read("q1_q17_dr25_koi.tab_selectcols_new.csv", header=19, nullable=true)
Q1Q17_DR25_stellar = CSV.read("q1_q17_dr25_stellar_koi.tab_selectcols.csv", header=30, nullable=true)
Q1Q17_DR25_stellar_all = CSV.read("q1_q17_dr25_stellar_koi.tab_all.csv", header=1, nullable=true)

N_Kepler_targets = sum((Q1Q17_DR25_stellar_all[:teff] .> 4000.) .& (Q1Q17_DR25_stellar_all[:teff] .< 7000.) .& (Q1Q17_DR25_stellar_all[:logg] .> 4.))
println("Total number of Kepler targets satisfying our cuts: ", N_Kepler_targets)

table_confirmed = Q1Q17_DR25[(Q1Q17_DR25[:koi_disposition] .== "CONFIRMED") .| (Q1Q17_DR25[:koi_disposition] .== "CANDIDATE"), :] #Table containing only the confirmed and candidate objects
table_stellar = Q1Q17_DR25_stellar

table_confirmed = table_confirmed[(table_confirmed[:koi_period] .> 5.) .& (table_confirmed[:koi_period] .< 300.), :] #to make additional cuts in period P to be comparable to our simulated sample
table_confirmed = table_confirmed[(table_confirmed[:koi_prad] .> 0.5) .& (table_confirmed[:koi_prad] .< 10.) .& (.~ismissing.(table_confirmed[:koi_prad])), :] #to make additional cuts in planetary radii to be comparable to our simulated sample

#To make cuts based on stellar properties of T_eff and logg:
teff_confirmed = zeros(Int64, size(table_confirmed,1)) #list to be filled with the T_eff (K) for the objects
logg_confirmed = zeros(Float64, size(table_confirmed,1)) #list to be filled with the logg(cgs) for the objects
cdpp5_confirmed = zeros(Float64, size(table_confirmed,1)) #list to be filled with the RMS CDPP 5h values for the objects
cdpp6_confirmed = zeros(Float64, size(table_confirmed,1)) #list to be filled with the RMS CDPP 6h values for the objects
for (i,KepID) in enumerate(table_confirmed[:kepid])
    teff_confirmed[i] = table_stellar[:teff][table_stellar[:kepid] .== KepID, :][1]
    logg_confirmed[i] = table_stellar[:logg][table_stellar[:kepid] .== KepID, :][1]
    cdpp5_confirmed[i] = table_stellar[:rrmscdpp05p0][table_stellar[:kepid] .== KepID, :][1]
    cdpp6_confirmed[i] = table_stellar[:rrmscdpp06p0][table_stellar[:kepid] .== KepID, :][1]
end

cdpp_cut = 250.
println("Fraction of CONFIRMED and CANDIDATE planets after CDPP cut: ", size(table_confirmed[(teff_confirmed .> 4000.) .& (teff_confirmed .< 7000.) .& (logg_confirmed .> 4.) .& (cdpp5_confirmed .< cdpp_cut), :], 1), "/", size(table_confirmed[(teff_confirmed .> 4000.) .& (teff_confirmed .< 7000.) .& (logg_confirmed .> 4.), :], 1))

table_confirmed = table_confirmed[(teff_confirmed .> 4000.) .& (teff_confirmed .< 7000.) .& (logg_confirmed .> 4.) .& (cdpp5_confirmed .< cdpp_cut), :]


KOI_systems = [x[1:6] for x in table_confirmed[:kepoi_name]]
checked_bools = zeros(size(table_confirmed,1)) #0's denote KOI that were not checked yet; 1's denote already checked KOI

M_confirmed = Int64[] #list to be filled with the planet multiplicities of the systems
R_confirmed = Float64[] #list to be filled with period ratios of adjacent planet pairs
xi_confirmed = Float64[] #list to be filled with the period-normalized transit duration ratios of adjacent planet pairs
D_ratio_confirmed = Float64[] #list to be filled with the transit depth ratios of adjacent planet pairs
P_confirmed = table_confirmed[:koi_period] #array of the periods (days)
P_confirmed = collect(skipmissing(P_confirmed))
t_D_confirmed = table_confirmed[:koi_duration] #array of the transit durations (hrs)
t_D_confirmed = collect(skipmissing(t_D_confirmed))
D_confirmed = table_confirmed[:koi_depth]/(1e6) #array of the transit depths (fraction)
D_confirmed = D_confirmed[ismissing.(D_confirmed) .== false] #to get rid of NA values

for i in 1:length(KOI_systems)
    if checked_bools[i] == 0 #if the KOI has not been checked (included while looking at another planet in the same system)
        system_i = (1:length(KOI_systems))[KOI_systems .== KOI_systems[i]]
        checked_bools[system_i] = 1

        #To get the periods and transit durations in this system:
        system_P = table_confirmed[:koi_period][system_i] #periods of all the planets in this system
        system_t_D = table_confirmed[:koi_duration][system_i] #transit durations of all the planets in this system
        system_D = table_confirmed[:koi_depth][system_i] #transit depths (in ppm) of all the planets in this system
        system_sort_i = sortperm(system_P) #indices that would sort the periods of the planets in this system
        system_P = system_P[system_sort_i] #periods of all the planets in this system, sorted
        system_t_D = system_t_D[system_sort_i] #transit durations of all the planets in this system, sorted by period
        system_D = system_D[system_sort_i] #transit depths of all the planets in this system, sorted by period

        #To count the total number of planets in this system:
        push!(M_confirmed, length(system_P))

        #To compute the period ratios, period-normalized transit duration ratios, and transit depth ratios in this system:
        system_R = system_P[2:end]./system_P[1:end-1] #period ratios of all the adjacent planet pairs in this system
        system_D_ratio = system_D[2:end]./system_D[1:end-1] #transit depth ratios of all the adjacent planet pairs in this system
        system_xi = (system_t_D[1:end-1]./system_t_D[2:end]).*(system_P[2:end]./system_P[1:end-1]).^(1//3) #period-normalized transit duration ratios of all the adjacent planet pairs in this system

        append!(R_confirmed, system_R)
        append!(D_ratio_confirmed, system_D_ratio)
        append!(xi_confirmed, system_xi)
    end
end





##### To define functions for calculating the KS distances:

function calc_distance(ss1::ExoplanetsSysSim.CatalogSummaryStatistics, ss2::ExoplanetsSysSim.CatalogSummaryStatistics)
    #This function calculates the total KS distance between two populations generated by our model

    M_cat_obs1 = ones(Int64,0) #array to be filled with the number of transiting planets in each simulated system for ss1
    M_cat_obs2 = ones(Int64,0) #array to be filled with the number of transiting planets in each simulated system for ss2
    for k in 1:get_int(sim_param,"max_tranets_in_sys")
        append!(M_cat_obs1, k*ones(Int64, ss1.stat["num n-tranet systems"][k]))
        append!(M_cat_obs2, k*ones(Int64, ss2.stat["num n-tranet systems"][k]))
    end

    d = Array{Float64}(8)
    d[1] = abs(ss1.stat["num_tranets"]/ss1.stat["num targets"] - ss2.stat["num_tranets"]/ss2.stat["num targets"])
    d[2] = ksstats_ints(M_cat_obs1, M_cat_obs2)[5]
    d[3] = ksstats(ss1.stat["P list"], ss2.stat["P list"])[5]
    d[4] = ksstats(ss1.stat["period_ratio_list"], ss2.stat["period_ratio_list"])[5]
    d[5] = ksstats(ss1.stat["duration list"], ss2.stat["duration list"])[5]
    d[6] = ksstats(ss1.stat["duration_ratio_list"], ss2.stat["duration_ratio_list"])[5]
    d[7] = ksstats(ss1.stat["depth list"], ss2.stat["depth list"])[5]
    d[8] = ksstats(ss1.stat["radius_ratio_list"], ss2.stat["radius_ratio_list"])[5]

    println("Distances: ", d, [sum(d)])
    println(f, "Dist: ", d, [sum(d)]) #to write the distances to file
    return sum(d)
end

function calc_distance_Kepler(ss1::ExoplanetsSysSim.CatalogSummaryStatistics)
    #This function calculates the total KS distance between a population generated by our model and the actual Kepler population (which must be loaded in)

    M_cat_obs = ones(Int64,0) #array to be filled with the number of transiting planets in each simulated system
    for k in 1:get_int(sim_param,"max_tranets_in_sys")
        append!(M_cat_obs, k*ones(Int64, ss1.stat["num n-tranet systems"][k]))
    end

    d = Array{Float64}(8)
    d[1] = abs(ss1.stat["num_tranets"]/ss1.stat["num targets"] - length(P_confirmed)/N_Kepler_targets)
    d[2] = ksstats_ints(M_cat_obs, M_confirmed)[5]
    d[3] = ksstats(ss1.stat["P list"], P_confirmed)[5]
    d[4] = ksstats(ss1.stat["period_ratio_list"], R_confirmed)[5]
    d[5] = ksstats(ss1.stat["duration list"].*24, t_D_confirmed)[5] #transit durations in simulations are in days, while in the Kepler catalog are in hours
    d[6] = ksstats(ss1.stat["duration_ratio_list"], xi_confirmed)[5]
    d[7] = ksstats(ss1.stat["depth list"], D_confirmed)[5]
    d[8] = ksstats(ss1.stat["radius_ratio_list"].^2, D_ratio_confirmed)[5] #simulations save radius ratios while we computed transit duration ratios from the Kepler catalog

    println("Distances: ", d, [sum(d)])
    println(f, "Dist: ", d, [sum(d)]) #to write the distances to file
    return sum(d)
end


function target_function(active_param::Vector)
    #This function takes in the values of the active model parameters and outputs the total KS distance compared to a reference simulated catalog

    global sim_param, summary_stat_ref
    sim_param_here = deepcopy(sim_param)
    ExoplanetsSysSim.update_sim_param_from_vector!(active_param,sim_param_here)
    cat_phys = generate_kepler_physical_catalog(sim_param_here)
    cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param_here)
    cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param_here)
    summary_stat = calc_summary_stats_model(cat_obs,sim_param_here)

    dist = calc_distance(summary_stat,summary_stat_ref)
end

function target_function_Kepler(active_param::Vector)
    #This function takes in the values of the active model parameters and outputs the total KS distance compared to the actual Kepler population

    global sim_param
    println("Active parameter values: ", active_param)
    println(f, "Active_params: ", active_param) #to write the params to file

    sim_param_here = deepcopy(sim_param)
    ExoplanetsSysSim.update_sim_param_from_vector!(active_param,sim_param_here)
    cat_phys = generate_kepler_physical_catalog(sim_param_here)
    cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param_here)
    cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param_here)
    summary_stat = calc_summary_stats_model(cat_obs,sim_param_here)

    dist = calc_distance_Kepler(summary_stat)
end





##### To draw the initial values of the active parameters randomly within a search range:

active_param_keys = ["break_radius", "log_rate_clusters", "log_rate_planets_per_cluster", "mr_power_index", "num_mutual_hill_radii", "power_law_P", "power_law_r1", "power_law_r2", "sigma_hk", "sigma_incl", "sigma_incl_near_mmr", "sigma_log_radius_in_cluster", "sigma_logperiod_per_pl_in_cluster"]
active_params_box = [(0.5*ExoplanetsSysSim.earth_radius, 10.*ExoplanetsSysSim.earth_radius), (log(1.), log(3.)), (log(1.), log(3.)), (1., 4.), (5., 20.), (-1., 1.), (-5., 3.), (-5., 3.), (0., 0.2), (0., 5.), (0., 5.), (0., 0.5), (0., 0.5)] #search ranges for all of the active parameters

#To randomly draw (uniformly) a value for each active model parameter within its search range:
for (i,param_key) in enumerate(active_param_keys)
    active_param_draw = active_params_box[i][1] + (active_params_box[i][2] - active_params_box[i][1])*rand(1)
    add_param_active(sim_param,param_key,active_param_draw[1])
end





##### To start saving the model iterations in the optimization into a file:

model_name = "Clustered_P_R_broken_R_optimization"
optimization_number = "" #if want to run on the cluster with random initial active parameters: "_random"*ARGS[1]
max_evals = 5
file_name = model_name*optimization_number*"_targs"*string(get_int(sim_param,"num_targets_sim_pass_one"))*"_evals"*string(max_evals)*".txt"

f = open(file_name, "w")
println(f, "# All initial parameters:")
write_model_params(f, sim_param)





##### To run the same model multiple times to see how it compares to the simulated catalog:

tic()
println("# Active parameters: ", make_vector_of_active_param_keys(sim_param))
println(f, "# Active parameters: ", make_vector_of_active_param_keys(sim_param))
active_param_true = make_vector_of_sim_param(sim_param)
println("# True values: ", active_param_true)
println(f, "# Starting active parameter values: ", active_param_true)
println(f, "# Format: Dist: [distances][total distance]")
num_eval = 5 #20
results = map(x->target_function(active_param_true), 1:num_eval)
mean_dist = mean(results)
rms_dist = std(results)
println("# Distance using true values: ", mean_dist, " +/- ",rms_dist)
println(f, "# Distance using active parameter values: ", mean_dist, " +/- ", rms_dist)
t_elapsed = toc()

println(f, "# elapsed time: ", t_elapsed, " seconds")
println(f, "#")






##### To run the same model multiple times to see how it compares to the actual Kepler catalog:
#=
tic()
println("# Active parameters: ", make_vector_of_active_param_keys(sim_param))
active_param_true = make_vector_of_sim_param(sim_param)
println("# Active parameter values: ", active_param_true)
num_eval = 10
results = map(x->target_function_Kepler(active_param_true), 1:num_eval)
mean_dist = mean(results)
rms_dist = std(results)
println("# Distance using default values: ", mean_dist, " +/- ",rms_dist)
toc()
=#




##### To use automated optimization routines to optimize the active model parameters:

active_param = 2*active_param_true
#=
using Optim   # I didn't have such good results with my initial attempts with BFGS, so I moved on
opt_result = optimize(target_function, active_param, method=BFGS(), f_tol=rms_dist, allow_f_increases=true, show_trace=true)
=#

# Pkg.add("BlackBoxOptim")       # only need to do these once
# Pkg.checkout("BlackBoxOptim")  # needed to get the lastest version
using BlackBoxOptim              # see https://github.com/robertfeldt/BlackBoxOptim.jl for documentation

# Eventually, we should make these physically motivated limits when we don't know the true values
#sr = [(active_param_true[i] > 0)? (active_param_true[i]/2,active_param_true[i]*2): (active_param_true[i]*2,active_param_true[i]/2) for i in 1:length(active_param)]

println(f, "# Optimization active parameters search bounds: ", active_params_box)
println(f, "# Format: Active_params: [active parameter values]")
println(f, "# Format: Dist: [distances][total distance]")

# May want to experiment with different algorithms, number of evaluations, population size, etc.
tic()
#opt_result = bboptimize(target_function_Kepler; SearchRange = active_params_box, NumDimensions = length(active_param), Method = :adaptive_de_rand_1_bin_radiuslimited, PopulationSize = length(active_param)*4, MaxFuncEvals = 20, TraceMode = :verbose)

opt_result = bboptimize(target_function_Kepler; SearchRange = active_params_box, NumDimensions = length(active_param), Method = :adaptive_de_rand_1_bin_radiuslimited, PopulationSize = length(active_param)*4, MaxFuncEvals = max_evals, TargetFitness = mean_dist, FitnessTolerance = rms_dist, TraceMode = :verbose)
t_elapsed = toc()

println(f, "# best_candidate: ", best_candidate(opt_result))
println(f, "# best_fitness: ", best_fitness(opt_result))
println(f, "# elapsed time: ", t_elapsed, " seconds")
close(f)


