import DataArrays.skipmissing

include("clusters.jl")

sim_param = setup_sim_param_model()
add_param_fixed(sim_param,"num_targets_sim_pass_one",150060) #9)   # For "observed" data, use a realistic number of targets (after any cuts you want to perform)
add_param_fixed(sim_param,"max_incl_sys",80.0) #degrees; 0 (deg) for isotropic system inclinations; set closer to 90 (deg) for more transiting systems

const max_incl_sys = get_real(sim_param,"max_incl_sys")
cos_factor = cos(max_incl_sys*pi/180) #factor to divide the number of targets in simulation by to get the actual number of targets needed (with an isotropic distribution of system inclinations) to produce as many transiting systems for a single observer





##### To load the DR25 catalog and compute arrays of multiplicities, periods, period ratios, transit durations, transit depths, period-normalized transit duration ratios (xi), and transit depth ratios:

Q1Q17_DR25 = CSV.read("q1_q17_dr25_koi.tab_selectcols_new.csv", header=19, nullable=true)
Q1Q17_DR25_stellar = CSV.read("q1_q17_dr25_stellar_koi.tab_selectcols.csv", header=30, nullable=true)
Q1Q17_DR25_stellar_all = CSV.read("q1_q17_dr25_stellar_koi.tab_all.csv", header=1, nullable=true)

N_Kepler_targets = sum((Q1Q17_DR25_stellar_all[:teff] .> 4000.) .& (Q1Q17_DR25_stellar_all[:teff] .< 7000.) .& (Q1Q17_DR25_stellar_all[:logg] .> 4.))
println("Total number of Kepler targets satisfying our cuts: ", N_Kepler_targets)

table_confirmed = Q1Q17_DR25[(Q1Q17_DR25[:koi_disposition] .== "CONFIRMED") .| (Q1Q17_DR25[:koi_disposition] .== "CANDIDATE"), :] #Table containing only the confirmed and candidate objects
table_stellar = Q1Q17_DR25_stellar

table_confirmed = table_confirmed[(table_confirmed[:koi_period] .> get_real(sim_param,"min_period")) .& (table_confirmed[:koi_period] .< get_real(sim_param,"max_period")), :] #to make additional cuts in period P to be comparable to our simulated sample
table_confirmed = table_confirmed[(table_confirmed[:koi_prad] .> get_real(sim_param,"min_radius")/ExoplanetsSysSim.earth_radius) .& (table_confirmed[:koi_prad] .< get_real(sim_param,"max_radius")/ExoplanetsSysSim.earth_radius) .& (.~ismissing.(table_confirmed[:koi_prad])), :] #to make additional cuts in planetary radii to be comparable to our simulated sample

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





##### To define functions for calculating the distances:

function calc_distance(ss1::ExoplanetsSysSim.CatalogSummaryStatistics, ss2::ExoplanetsSysSim.CatalogSummaryStatistics, all_dist::Bool=false, save_dist::Bool=true)
    #This function calculates the total KS distance between two summary statistics (simulated observed catalogs).
    #If 'all_dist=true', the function outputs the individual distances in the distance function.
    #If 'save_dist=true', the function also saves the distances (individual and total) to a file (assuming file 'f' is open for writing).

    global cos_factor

    M_cat_obs1 = ones(Int64,0) #array to be filled with the number of transiting planets in each simulated system for ss1
    M_cat_obs2 = ones(Int64,0) #array to be filled with the number of transiting planets in each simulated system for ss2
    for k in 1:get_int(sim_param,"max_tranets_in_sys")
        append!(M_cat_obs1, k*ones(Int64, ss1.stat["num n-tranet systems"][k]))
        append!(M_cat_obs2, k*ones(Int64, ss2.stat["num n-tranet systems"][k]))
    end

    #To handle empty arrays:
    if sum(M_cat_obs1 .>= 2) < 2 | sum(M_cat_obs2 .>= 2) < 2 #need at least 2 multi-systems per catalog in order to be able to compute AD distances for distributions of ratios of observables
        println("Not enough observed multi-planet systems in one of the catalogs to compute the AD distance.")
        d = ones(Int64,8)*1e6

        println("Distances: ", d, [sum(d)])
        if save_dist
            println(f, "Dist_KS: ", d, [sum(d)])
            #println(f, "Dist_AD: ", d, [sum(d)])
        end

        if all_dist
            return d
        else
            return sum(d)
        end
    end

    #To compute the KS distances:
    d_KS = Array{Float64}(8)
    d_KS[1] = abs(ss1.stat["num_tranets"]/(ss1.stat["num targets"]/cos_factor) - ss2.stat["num_tranets"]/(ss2.stat["num targets"]/cos_factor))
    d_KS[2] = ksstats_ints(M_cat_obs1, M_cat_obs2)[5]
    d_KS[3] = ksstats(ss1.stat["P list"], ss2.stat["P list"])[5]
    d_KS[4] = ksstats(ss1.stat["period_ratio_list"], ss2.stat["period_ratio_list"])[5]
    d_KS[5] = ksstats(ss1.stat["duration list"], ss2.stat["duration list"])[5]
    d_KS[6] = ksstats(ss1.stat["duration_ratio_list"], ss2.stat["duration_ratio_list"])[5]
    d_KS[7] = ksstats(ss1.stat["depth list"], ss2.stat["depth list"])[5]
    d_KS[8] = ksstats(ss1.stat["radius_ratio_list"], ss2.stat["radius_ratio_list"])[5]

    #=
    #To compute the AD distances:
    d_AD = Array{Float64}(8)
    d_AD[1] = abs(ss1.stat["num_tranets"]/(ss1.stat["num targets"]/cos_factor) - ss2.stat["num_tranets"]/(ss2.stat["num targets"]/cos_factor))
    d_AD[2] = ksstats_ints(M_cat_obs1, M_cat_obs2)[5]
    d_AD[3] = ADstats(ss1.stat["P list"], ss2.stat["P list"])
    d_AD[4] = ADstats(ss1.stat["period_ratio_list"], ss2.stat["period_ratio_list"])
    d_AD[5] = ADstats(ss1.stat["duration list"], ss2.stat["duration list"])
    d_AD[6] = ADstats(ss1.stat["duration_ratio_list"], ss2.stat["duration_ratio_list"])
    d_AD[7] = ADstats(ss1.stat["depth list"], ss2.stat["depth list"])
    d_AD[8] = ADstats(ss1.stat["radius_ratio_list"], ss2.stat["radius_ratio_list"])
    =#

    #To print and/or write the distances to file:
    println("KS Distances: ", d_KS, [sum(d_KS)])
    #println("AD Distances: ", d_AD, [sum(d_AD)])
    if save_dist
        println(f, "Dist_KS: ", d_KS, [sum(d_KS)])
        #println(f, "Dist_AD: ", d_AD, [sum(d_AD)])
    end

    #To return the distances or total distance:
    if all_dist
        return d_KS
        #return d_AD
    else
        return sum(d_KS)
        #return sum(d_AD)
    end
end

function calc_distance_Kepler(ss1::ExoplanetsSysSim.CatalogSummaryStatistics, all_dist::Bool=false, save_dist::Bool=true)
    #This function calculates the total KS distance between a population generated by our model and the actual Kepler population (which must be loaded in).
    #If 'all_dist=true', the function outputs the individual distances in the distance function.
    #If 'save_dist=true', the function also saves the distances (individual and total) to a file (assuming file 'f' is open for writing).

    global cos_factor

    M_cat_obs = ones(Int64,0) #array to be filled with the number of transiting planets in each simulated system
    for k in 1:get_int(sim_param,"max_tranets_in_sys")
        append!(M_cat_obs, k*ones(Int64, ss1.stat["num n-tranet systems"][k]))
    end

    if sum(M_cat_obs .>= 2) < 2 #need at least 2 observed multi-systems in order to be able to compute AD distances for distributions of ratios of observables
        println("Not enough observed multi-planet systems in the simulated catalog.")
        d = ones(Int64,8)*1e6

        println("Distances: ", d, [sum(d)])
        if save_dist
            println(f, "Dist_KS: ", d, [sum(d)])
            #println(f, "Dist_AD: ", d, [sum(d)])
        end

        if all_dist
            return d
        else
            return sum(d)
        end
    end

    #To compute the KS distances:
    d_KS = Array{Float64}(8)
    d_KS[1] = abs(ss1.stat["num_tranets"]/(ss1.stat["num targets"]/cos_factor) - length(P_confirmed)/N_Kepler_targets)
    d_KS[2] = ksstats_ints(M_cat_obs, M_confirmed)[5]
    d_KS[3] = ksstats(ss1.stat["P list"], P_confirmed)[5]
    d_KS[4] = ksstats(ss1.stat["period_ratio_list"], R_confirmed)[5]
    d_KS[5] = ksstats(ss1.stat["duration list"].*24, t_D_confirmed)[5] #transit durations in simulations are in days, while in the Kepler catalog are in hours
    d_KS[6] = ksstats(ss1.stat["duration_ratio_list"], xi_confirmed)[5]
    d_KS[7] = ksstats(ss1.stat["depth list"], D_confirmed)[5]
    d_KS[8] = ksstats(ss1.stat["radius_ratio_list"].^2, D_ratio_confirmed)[5] #simulations save radius ratios while we computed transit duration ratios from the Kepler catalog

    #=
    #To compute the AD distances:
    d_AD = Array{Float64}(8)
    d_AD[1] = abs(ss1.stat["num_tranets"]/(ss1.stat["num targets"]/cos_factor) - length(P_confirmed)/N_Kepler_targets)
    d_AD[2] = ksstats_ints(M_cat_obs, M_confirmed)[5]
    d_AD[3] = ADstats(ss1.stat["P list"], P_confirmed)
    d_AD[4] = ADstats(ss1.stat["period_ratio_list"], R_confirmed)
    d_AD[5] = ADstats(ss1.stat["duration list"].*24, t_D_confirmed) #transit durations in simulations are in days, while in the Kepler catalog are in hours
    d_AD[6] = ADstats(ss1.stat["duration_ratio_list"], xi_confirmed)
    d_AD[7] = ADstats(ss1.stat["depth list"], D_confirmed)
    d_AD[8] = ADstats(ss1.stat["radius_ratio_list"].^2, D_ratio_confirmed) #simulations save radius ratios while we computed transit duration ratios from the Kepler catalog
    =#

    #To print and/or write the distances to file:
    println("KS Distances: ", d_KS, [sum(d_KS)])
    #println("AD Distances: ", d_AD, [sum(d_AD)])
    if save_dist
        println(f, "Dist_KS: ", d_KS, [sum(d_KS)])
        #println(f, "Dist_AD: ", d_AD, [sum(d_AD)])
    end

    #To return the distances or total distance:
    if all_dist
        return d_KS
        #return d_AD
    else
        return sum(d_KS)
        #return sum(d_AD)
    end
end


function target_function(active_param::Vector, all_dist::Bool=false, save_dist::Bool=true)
    #This function takes in the values of the active model parameters, generates a simulated observed catalog, and computes the distance function compared to a reference simulated catalog.
    #If 'all_dist=true', the function outputs the individual distances in the distance function.
    #If 'save_dist=true', the function also saves the distances (individual and total) to a file (assuming file 'f' is open for writing).

    println("Active parameter values:", active_param)
    println(f, "Active_params: ", active_param) #if we also want to write the params to file

    global sim_param, summary_stat_ref
    sim_param_here = deepcopy(sim_param)
    ExoplanetsSysSim.update_sim_param_from_vector!(active_param,sim_param_here)
    cat_phys = generate_kepler_physical_catalog(sim_param_here)
    cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param_here)
    cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param_here)
    summary_stat = calc_summary_stats_model(cat_obs,sim_param_here)

    dist = calc_distance(summary_stat, summary_stat_ref, all_dist, save_dist)
    return dist
end

function target_function_Kepler(active_param::Vector, all_dist::Bool=false, save_dist::Bool=true)
    #This function takes in the values of the active model parameters, generates a simulated observed catalog, and computes the distance function compared to the actual Kepler population.
    #If 'all_dist=true', the function outputs the individual distances in the distance function.
    #If 'save_dist=true', the function also saves the distances (individual and total) to a file (assuming file 'f' is open for writing).

    println("Active parameter values: ", active_param)
    if save_dist
        println(f, "Active_params: ", active_param) #to write the params to file
    end

    global sim_param
    sim_param_here = deepcopy(sim_param)
    ExoplanetsSysSim.update_sim_param_from_vector!(active_param,sim_param_here)
    cat_phys = generate_kepler_physical_catalog(sim_param_here)
    cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param_here)
    cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param_here)
    summary_stat = calc_summary_stats_model(cat_obs,sim_param_here)

    dist = calc_distance_Kepler(summary_stat, all_dist, save_dist)
    return dist
end

function target_function_weighted(active_param::Vector, all_dist::Bool=false, save_dist::Bool=true)
    #This function takes in the values of the active model parameters, generates a simulated observed catalog, and computes the weighted distance function (assuming an array of weights 'weights' is a global variable) compared to a reference simulated catalog.
    #If 'all_dist=true', the function outputs the individual distances (weighted) in the distance function.
    #If 'save_dist=true', the function also saves the distances (unweighted and weighted, individual and total) to a file (assuming file 'f' is open for writing).

    global weights

    dist = target_function(active_param, true, save_dist)
    weighted_dist = dist./weights
    #used_dist = weighted_dist[1:6] #choose a subset of the distances to pass into the optimizer

    println("Weighted distances: ", weighted_dist, [sum(weighted_dist)])
    if save_dist
        println(f, "Dist_weighted: ", weighted_dist, [sum(weighted_dist)])
    end

    if all_dist
        return weighted_dist
    else
        return sum(weighted_dist)
        #return sum(used_dist)
    end
end

function target_function_Kepler_weighted(active_param::Vector, all_dist::Bool=false, save_dist::Bool=true)
    #This function takes in the values of the active model parameters, generates a simulated observed catalog, and computes the weighted distance function (assuming an array of weights 'weights' is a global variable) compared to the actual Kepler population.
    #If 'all_dist=true', the function outputs the individual distances (weighted) in the distance function.
    #If 'save_dist=true', the function also saves the distances (unweighted and weighted, individual and total) to a file (assuming file 'f' is open for writing).

    global weights

    dist = target_function_Kepler(active_param, true, save_dist)
    weighted_dist = dist./weights
    #used_dist = weighted_dist[1:6] #choose a subset of the distances to pass into the optimizer

    println("Weighted distances: ", weighted_dist, [sum(weighted_dist)])
    if save_dist
        println(f, "Dist_weighted: ", weighted_dist, [sum(weighted_dist)])
    end

    if all_dist
        return weighted_dist
    else
        return sum(weighted_dist)
        #return sum(used_dist)
    end
end





##### To start saving the model iterations in the optimization into a file:

model_name = "Clustered_P_R_broken_R_simulated_optimization_blocks"
optimization_number = "" #if want to run on the cluster with random initial active parameters: "_random"*ARGS[1]
max_evals = 100
cycles = 3
file_name = model_name*optimization_number*"_targs"*string(get_int(sim_param,"num_targets_sim_pass_one"))*"_evals"*string(max_evals)*"_cycles"*string(cycles)*".txt"

f = open(file_name, "w")
println(f, "# All initial parameters (for estimating stochastic noise in a perfect model):")
write_model_params(f, sim_param)





##### To run the same model multiple times to see how it compares to a simulated catalog with the same parameters:

#To generate a simulated catalog to fit to:
cat_phys = generate_kepler_physical_catalog(sim_param)
cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param)
cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param)
summary_stat_ref = calc_summary_stats_model(cat_obs,sim_param)

tic()
active_param_true = make_vector_of_sim_param(sim_param)
println("# True values: ", active_param_true)
println(f, "# Format: Dist: [distances][total distance]")

num_eval = 20 #20
dists_true = zeros(num_eval,8)
for i in 1:num_eval
    dists_true[i,:] = target_function(active_param_true, true)
end
mean_dists = transpose(mean(dists_true,1))[:,] #array of mean distances for each individual distance
mean_dist = mean(sum(dists_true, 2)) #mean total distance
#std_dists = transpose(std(dists_true, 1))[:,] #array of std distances for each individual distance
rms_dists = transpose(sqrt.(mean(dists_true.^2, 1)))[:,] #array of rms (std around 0)  distances for each individual distance
std_dist = std(sum(dists_true, 2)) #std of total distance
rms_dist = sqrt(mean(sum(dists_true, 2).^2)) #rms (std around 0) of total distance; should be similar to the mean total distance

weights = rms_dists #to use the array 'rms_dists' as the weights for the individual distances
weighted_dists_true = zeros(num_eval,8)
for i in 1:num_eval
    weighted_dists_true[i,:] = dists_true[i,:]./rms_dists
end
mean_weighted_dists = transpose(mean(weighted_dists_true,1))[:,] #array of mean weighted distances for each individual distance
mean_weighted_dist = mean(sum(weighted_dists_true,2)) #mean weighted total distance
std_weighted_dist = std(sum(weighted_dists_true,2)) #std of weighted total distance

println("Mean dists: ", mean_dists)
println("Rms dists: ", rms_dists)
println("Mean weighted dists: ", mean_weighted_dists)
println("Distance using true values: ", mean_dist, " +/- ", std_dist)
println("Weighted distance using true values: ", mean_weighted_dist, " +/- ", std_weighted_dist)
println(f, "Mean: ", mean_dists, [mean_dist])
println(f, "Rms: ", rms_dists, [rms_dist])
println(f, "Mean weighted dists: ", mean_weighted_dists, [mean_weighted_dist])
println(f, "# Distance using true values (default parameter values): ", mean_dist, " +/- ", std_dist)
println(f, "# Weighted distance using true values (default parameter values): ", mean_weighted_dist, " +/- ", std_weighted_dist)
t_elapsed = toc()
println(f, "# elapsed time: ", t_elapsed, " seconds")
println(f, "#")





##### To set up a block-optimization routine:

#To specify the blocks of random parameters, their search ranges, and to draw initial values:
active_param_keys = [["log_rate_clusters", "log_rate_planets_per_cluster", "num_mutual_hill_radii", "power_law_P", "sigma_logperiod_per_pl_in_cluster"], ["break_radius", "mr_power_index", "power_law_r1", "power_law_r2", "sigma_log_radius_in_cluster"], ["sigma_hk", "sigma_incl", "sigma_incl_near_mmr"]]
    #[["log_rate_clusters", "log_rate_planets_per_cluster", "num_mutual_hill_radii", "power_law_P", "sigma_logperiod_per_pl_in_cluster"], ["break_radius", "mr_power_index", "power_law_r1", "power_law_r2", "sigma_log_radius_in_cluster"], ["sigma_hk", "sigma_incl", "sigma_incl_near_mmr"]]
active_params_box = [[(log(1.), log(5.)), (log(1.), log(5.)), (3., 20.), (-0.5, 1.5), (0., 0.3)], [(0.5*ExoplanetsSysSim.earth_radius, 10.*ExoplanetsSysSim.earth_radius),(1., 4.), (-6., 0.), (-6., 0.), (0.1, 1.0)], [(0., 0.1), (0., 5.), (0., 5.)]] #search ranges for all of the active parameters
    #[[(log(1.), log(5.)), (log(1.), log(5.)), (3., 20.), (-0.5, 1.5), (0., 0.3)], [(0.5*ExoplanetsSysSim.earth_radius, 10.*ExoplanetsSysSim.earth_radius),(1., 4.), (-6., 0.), (-6., 0.), (0.1, 1.0)], [(0., 0.1), (0., 5.), (0., 5.)]] #search ranges for all of the active parameters

active_param_draws = [Float64[], Float64[], Float64[]] #[Float64[], Float64[], Float64[]]
for (i,active_param_block) in enumerate(active_param_keys)
    for (j,param_key) in enumerate(active_param_block)
        active_param_draw = active_params_box[i][j][1] + (active_params_box[i][j][2] - active_params_box[i][j][1])*rand(1)
        append!(active_param_draws[i], active_param_draw)
        add_param_fixed(sim_param, param_key, active_param_draw[1])
    end
end

println("# All active parameters: ", active_param_keys)
println(f, "# All active parameters: ", active_param_keys)
println(f, "# Starting active parameter values: ", active_param_draws)
println(f, "# Optimization active parameters search bounds: ", active_params_box)
println(f, "# Method: adaptive_de_rand_1_bin_radiuslimited")
println(f, "# Format: Active_params: [active parameter values]")
println(f, "# Format: Dist: [distances][total distance]")
println(f, "# Format: Dist_weighted: [weighted distances][total weighted distance]")
println(f, "# Distances used: Dist_KS (all, weighted by rms)") #Edit this line to specify which distances were actually used in the optimizer!

#To start the block-optimization:

# Pkg.add("BlackBoxOptim")       # only need to do these once
# Pkg.checkout("BlackBoxOptim")  # needed to get the lastest version
using BlackBoxOptim              # see https://github.com/robertfeldt/BlackBoxOptim.jl for documentation

tic()
for c in 1:cycles
    println("# Cycle: ", c)
    println(f, "#")
    println(f, "# Cycle: ", c)
    println(f, "#")

    tic()
    for (i,active_param_block) in enumerate(active_param_keys)
        println("# Active parameters: ", active_param_block)
        println(f, "# Active parameters: ", active_param_block)

        #To make the current block of parameters active:
        for (j,param_key) in enumerate(active_param_block)
            add_param_active(sim_param, param_key, active_param_draws[i][j])
        end

        tic()
        opt_result = bboptimize(target_function_weighted; SearchRange = active_params_box[i], NumDimensions = length(active_param_block), Method = :adaptive_de_rand_1_bin_radiuslimited, PopulationSize = length(active_param_true)*4, MaxFuncEvals = max_evals, TargetFitness = mean_weighted_dist, FitnessTolerance = std_weighted_dist, TraceMode = :verbose)
        t_block = toc()
        active_params_best = best_candidate(opt_result)
        best_distance = best_fitness(opt_result)

        println(f, "# Best parameters block: ", active_params_best)
        println(f, "# Best fitness block: ", best_distance)
        println(f, "# Block elapsed time: ", t_block, " seconds")
        println(f, "#")

        #To make the block of parameters fixed at the best values found:
        for (j,param_key) in enumerate(active_param_block)
            add_param_fixed(sim_param, param_key, active_params_best[j])
        end
    end
    t_cycle = toc()

    println(f, "# Cycle elapsed time: ", t_cycle, " seconds")
end
t_optim = toc()

println(f, "#")
println(f, "# Total elapsed time: ", t_optim, " seconds")
close(f)




