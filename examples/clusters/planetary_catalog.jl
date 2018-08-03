import DataArrays.skipmissing

sim_param = setup_sim_param_model()





##### To load the DR25 KOI catalog and apply the necessary cuts:

Q1Q17_DR25 = CSV.read("q1_q17_dr25_koi.tab_selectcols_new.csv", header=19, nullable=true)


#If using our own cuts applied to the stellar and koi catalogs:

#=
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
=#


#If using a stellar table with cuts already made and just matching kepids to get a koi catalog:

stellar_catalog = setup_star_table_christiansen(sim_param)
N_Kepler_targets = get_int(sim_param,"num_kepler_targets")

table_confirmed = Q1Q17_DR25[(Q1Q17_DR25[:koi_disposition] .== "CONFIRMED") .| (Q1Q17_DR25[:koi_disposition] .== "CANDIDATE"), :] #Table containing only the confirmed and candidate objects

table_confirmed = table_confirmed[(table_confirmed[:koi_period] .> get_real(sim_param,"min_period")) .& (table_confirmed[:koi_period] .< get_real(sim_param,"max_period")), :] #to make additional cuts in period P to be comparable to our simulated sample
table_confirmed = table_confirmed[(table_confirmed[:koi_prad] .> get_real(sim_param,"min_radius")/ExoplanetsSysSim.earth_radius) .& (table_confirmed[:koi_prad] .< get_real(sim_param,"max_radius")/ExoplanetsSysSim.earth_radius) .& (.~ismissing.(table_confirmed[:koi_prad])), :] #to make additional cuts in planetary radii to be comparable to our simulated sample

in_stellar_catalog = [] #will be filled with booleans indicating whether each koi in 'table_confirmed' is found in the 'stellar_catalog' or not
for i in 1:length(table_confirmed[:kepid])
    if any(x->x==table_confirmed[i,:kepid], stellar_catalog[:kepid])
        push!(in_stellar_catalog, true)
    else
        push!(in_stellar_catalog, false)
    end
end
in_stellar_catalog_indices = find(in_stellar_catalog)

#=
#If we want to write the kep_oi's of the remaining koi catalog to a file:
f = open("kepoi_names.txt", "w")
println(f, table_confirmed[in_stellar_catalog_indices, :kepoi_name])
close(f)
=#

table_confirmed = table_confirmed[in_stellar_catalog_indices, :]





#####To compute arrays of the observables (multiplicities, periods, period ratios, transit durations, transit depths, period-normalized transit duration ratios (xi), and transit depth ratios) from the remaining sample of planets:

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


