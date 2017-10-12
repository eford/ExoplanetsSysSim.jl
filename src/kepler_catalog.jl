## ExoplanetsSysSim/src/kepler_catalog.jl
## (c) 2015 Eric B. Ford

#using ExoplanetsSysSim
#using DataFrames
using JLD

#if VERSION >= v"0.5-"
#  import Compat: UTF8String, ASCIIString
#end

type KeplerPhysicalCatalog
  target::Array{KeplerTarget,1}
end
#KeplerPhysicalCatalog() = KeplerPhysicalCatalog([])

type KeplerObsCatalog
  target::Array{KeplerTargetObs,1}
end
#KeplerObsCatalog() = KeplerObsCatalog([])

function generate_kepler_physical_catalog(sim_param::SimParam)
   if haskey(sim_param,"stellar_catalog")
      star_tab_func = get_function(sim_param, "star_table_setup")
      star_tab_func(sim_param)
   end
   num_sys = get_int(sim_param,"num_targets_sim_pass_one")
   generate_kepler_target = get_function(sim_param,"generate_kepler_target")
   target_list = Array{KeplerTarget}(num_sys)
   map!(x->generate_kepler_target(sim_param), target_list, 1:num_sys )
   return KeplerPhysicalCatalog(target_list)
end

function observe_kepler_targets_sky_avg(input::KeplerPhysicalCatalog, sim_param::SimParam )
  calc_target_obs = get_function(sim_param,"calc_target_obs_sky_ave")
  return observe_kepler_targets(calc_target_obs, input, sim_param)
end

function observe_kepler_targets_single_obs(input::KeplerPhysicalCatalog, sim_param::SimParam )
  calc_target_obs = get_function(sim_param,"calc_target_obs_single_obs")
  return observe_kepler_targets(calc_target_obs, input, sim_param)
end

function observe_kepler_targets(calc_target_obs::Function, input::KeplerPhysicalCatalog, sim_param::SimParam )
  #calc_target_obs = get_function(sim_param,"calc_target_obs_sky_ave")
  #calc_target_obs = get_function(sim_param,"calc_target_obs_single_obs")
  output = KeplerObsCatalog([])
  if haskey(sim_param,"mem_kepler_target_obs")
     output.target = get(sim_param,"mem_kepler_target_obs",Array{KeplerTargetObs}(0) )
  end
  num_targets_sim_pass_one = get_int(sim_param,"num_targets_sim_pass_one")
  if length(output.target) < num_targets_sim_pass_one
     output.target = Array{KeplerTargetObs}(num_targets_sim_pass_one) 
  end
  #output.target = Array{KeplerTargetObs}(length(input.target) )  # Replaced to reduce memory allocation
  map!(x::KeplerTarget->calc_target_obs(x,sim_param)::KeplerTargetObs, output.target, input.target)
  resize!(output.target,length(input.target))
  return output
end

# Test if this planetary system has at least one planet that transits (assuming a single observer)
function select_targets_one_obs(ps::PlanetarySystemAbstract)
 for pl in 1:length(ps.orbit)
   ecc::Float64 = ps.orbit[pl].ecc
   incl::Float64 = ps.orbit[pl].incl
   a::Float64 = semimajor_axis(ps,pl)
   Rstar::Float64 = rsol_in_au*ps.star.radius
   if (Rstar > (a*(1-ecc)*(1+ecc))/(1+ecc*sin(ps.orbit[pl].omega))*cos(incl))
     return true
   end
 end
 return false
end

# Remove undetected planets from physical catalog
# TODO: OPT: Maybe create array of bools for which planets to keep, rather than splicing out non-detections?
function generate_obs_targets(cat_phys::KeplerPhysicalCatalog, sim_param::SimParam )
  for t in 1:length(cat_phys.target)
    for ps in 1:length(cat_phys.target[t].sys)
      kep_targ = cat_phys.target[t].sys[ps]
      for pl in length(kep_targ.orbit):-1:1    # Going in reverse since removing planets from end of list first is cheaper than starting at beginning
        ecc::Float64 = kep_targ.orbit[pl].ecc
	incl::Float64 = kep_targ.orbit[pl].incl
   	a::Float64 = semimajor_axis(kep_targ,pl)
   	Rstar::Float64 = rsol_in_au*kep_targ.star.radius
       
   	if (Rstar < (a*(1-ecc)*(1+ecc))/(1+ecc*sin(kep_targ.orbit[pl].omega))*cos(incl)) || (rand() > calc_prob_detect_if_transit(cat_phys.target[t], 1, pl, sim_param))
    	  splice!(cat_phys.target[t].sys[ps].orbit, pl)
	  splice!(cat_phys.target[t].sys[ps].planet, pl)
     	end
      end
    end
  end
  return cat_phys
end

function simulated_read_kepler_observations(sim_param::SimParam ) # TODO SCI:  IMPORTANT:  Eventually, replace this with a function to read data from input file (see koi_table.jl)
   if haskey(sim_param,"stellar_catalog")
      star_tab_func = get_function(sim_param, "star_table_setup")
      star_tab_func(sim_param)
   end
   num_sys = get_int(sim_param,"num_kepler_targets")
   generate_kepler_target = get_function(sim_param,"generate_kepler_target")
   target_list = Array{KeplerTarget}(num_sys)
   map!(x->generate_kepler_target(sim_param), target_list, 1:num_sys ) 

   cat_phys_cut = generate_obs_targets(KeplerPhysicalCatalog(target_list), sim_param)
   calc_target_obs = get_function(sim_param,"calc_target_obs_single_obs")
   output = KeplerObsCatalog([])
   output.target = map(x::KeplerTarget->calc_target_obs(x,sim_param)::KeplerTargetObs, cat_phys_cut.target)
   return output
end

# df_star is assumed to have fields kepid, mass and radius for all targets in the survey
function setup_actual_planet_candidate_catalog(df_star::DataFrame, sim_param::SimParam)
  local csv_data,kepid_idx,koi_period_idx,koi_time0bk_idx,koi_depth_idx,koi_duration_idx,koi_ror_idx,koi_period_err1_idx,koi_time0bk_err1_idx,koi_depth_err1_idx,koi_duration_err1_idx,koi_period_err2_idx,koi_time0bk_err2_idx,koi_depth_err2_idx,koi_duration_err2_idx 
  add_param_fixed(sim_param,"num_kepler_targets",num_usable_in_star_table())  # For "observed" catalog
  koi_catalog_file_in = convert(String,joinpath(Pkg.dir("ExoplanetsSysSim"), "data", convert(String,get(sim_param,"koi_catalog","q1_q17_dr25_koi.csv")) ) )

  if ismatch(r".jld$",koi_catalog_file_in)
  try 
    data = load(koi_catalog_file_in)
    csv_data = data["koi_catalog"]

    kepid_idx = 1
    koi_period_idx = 2
    koi_time0bk_idx = 3
    koi_depth_idx = 4
    koi_duration_idx = 5
    koi_ror_idx = 6
    koi_period_err1_idx = 7
    koi_time0bk_err1_idx = 8
    koi_depth_err1_idx = 9
    koi_duration_err1_idx = 10
    koi_period_err2_idx = 11
    koi_time0bk_err2_idx = 12
    koi_depth_err2_idx = 13
    koi_duration_err2_idx = 14
  catch
    error(string("# Failed to read koi catalog >",koi_catalog_file_in,"< in jld format."))
  end
  else
  try 
    (csv_data,csv_header) = readcsv(koi_catalog_file_in,header=true)
  
    # Lookup header columns, since DataFrames doesn't like this file
    kepid_idx = findfirst(x->x=="kepid",csv_header)
    koi_period_idx = findfirst(x->x=="koi_period",csv_header)
    koi_time0bk_idx = findfirst(x->x=="koi_time0bk",csv_header)
    koi_depth_idx = findfirst(x->x=="koi_depth",csv_header)
    koi_duration_idx = findfirst(x->x=="koi_duration",csv_header)
    koi_ror_idx = findfirst(x->x=="koi_ror",csv_header)
    koi_period_err1_idx = findfirst(x->x=="koi_period_err1",csv_header)
    koi_time0bk_err1_idx = findfirst(x->x=="koi_time0bk_err1",csv_header)
    koi_depth_err1_idx = findfirst(x->x=="koi_depth_err1",csv_header)
    koi_duration_err1_idx = findfirst(x->x=="koi_duration_err1",csv_header)
    koi_period_err2_idx = findfirst(x->x=="koi_period_err2",csv_header)
    koi_time0bk_err2_idx = findfirst(x->x=="koi_time0bk_err2",csv_header)
    koi_depth_err2_idx = findfirst(x->x=="koi_depth_err2",csv_header)
    koi_duration_err2_idx = findfirst(x->x=="koi_duration_err2",csv_header)
    koi_disposition_idx = findfirst(x->x=="koi_disposition",csv_header)
    koi_pdisposition_idx = findfirst(x->x=="koi_pdisposition",csv_header)
    # Choose which KOIs to keep
    #is_cand = (csv_data[:,koi_disposition_idx] .== "CONFIRMED") | (csv_data[:,koi_disposition_idx] .== "CANDIDATE")
    is_cand = (csv_data[:,koi_pdisposition_idx] .== "CANDIDATE")
  
    idx_keep = is_cand .& .!isna.(csv_data[:,koi_ror_idx]) .& ([typeof(x) for x in csv_data[:,koi_ror_idx]] .== Float64)
    idx_keep = idx_keep .& .!isna.(csv_data[:,koi_period_err1_idx]) .& ([typeof(x) for x in csv_data[:,koi_period_err1_idx]] .== Float64) # DR25 catalog missing uncertainties for some candidates
    csv_data = csv_data[idx_keep,:]
  catch
    error(string("# Failed to read koi catalog >",koi_catalog_file_in,"< in ascii format."))
  end
  end

  output = KeplerObsCatalog([])

  for (j,kepid) in enumerate(df_star[:kepid])
     plids = find(x->x==kepid,csv_data[:,kepid_idx])
	 num_pl = length(plids)
	 if num_pl ==0 continue end
	 target_obs = KeplerTargetObs(num_pl)
	 target_obs.star = ExoplanetsSysSim.StarObs(df_star[j,:radius],df_star[j,:mass])
	 for (plid,i) in enumerate(plids)
	   target_obs.obs[plid] = ExoplanetsSysSim.TransitPlanetObs(csv_data[i,koi_period_idx],csv_data[i,koi_time0bk_idx],csv_data[i,koi_depth_idx]/1.0e6,csv_data[i,koi_duration_idx])
           target_obs.sigma[plid] = ExoplanetsSysSim.TransitPlanetObs((abs(csv_data[i,koi_period_err1_idx])+abs(csv_data[i,koi_period_err2_idx]))/2,(abs(csv_data[i,koi_time0bk_err1_idx])+abs(csv_data[i,koi_time0bk_err2_idx]))/2,(abs(csv_data[i,koi_depth_err1_idx]/1.0e6)+abs(csv_data[i,koi_depth_err2_idx]/1.0e6))/2,(abs(csv_data[i,koi_duration_err1_idx])+abs(csv_data[i,koi_duration_err2_idx]))/2)
	   #target_obs.prob_detect = ExoplanetsSysSim.SimulatedSystemDetectionProbs{OneObserver}( ones(num_pl), ones(num_pl,num_pl), ones(num_pl), fill(Array{Int64}(0), 1) )  # Made line below to simplify calling
           target_obs.prob_detect = ExoplanetsSysSim.OneObserverSystemDetectionProbs(num_pl)
	 end	
      push!(output.target,target_obs)
  end
  return output
end

function setup_actual_planet_candidate_catalog_csv(df_star::DataFrame, sim_param::SimParam)
  add_param_fixed(sim_param,"num_kepler_targets",num_usable_in_star_table())  # For "observed" catalog
  koi_catalog_file_in = convert(String,joinpath(Pkg.dir("ExoplanetsSysSim"), "data", convert(String,get(sim_param,"koi_catalog","q1_q17_dr25_koi.csv")) ) )
  (csv_data,csv_header) = readcsv(koi_catalog_file_in,header=true)

  # Lookup header columns, since DataFrames doesn't like this file
  kepid_idx = findfirst(x->x=="kepid",csv_header)
  koi_period_idx = findfirst(x->x=="koi_period",csv_header)
  koi_time0bk_idx = findfirst(x->x=="koi_time0bk",csv_header)
  koi_depth_idx = findfirst(x->x=="koi_depth",csv_header)
  koi_duration_idx = findfirst(x->x=="koi_duration",csv_header)
  koi_ror_idx = findfirst(x->x=="koi_ror",csv_header)
  koi_period_err1_idx = findfirst(x->x=="koi_period_err1",csv_header)
  koi_time0bk_err1_idx = findfirst(x->x=="koi_time0bk_err1",csv_header)
  koi_depth_err1_idx = findfirst(x->x=="koi_depth_err1",csv_header)
  koi_duration_err1_idx = findfirst(x->x=="koi_duration_err1",csv_header)
  koi_period_err2_idx = findfirst(x->x=="koi_period_err2",csv_header)
  koi_time0bk_err2_idx = findfirst(x->x=="koi_time0bk_err2",csv_header)
  koi_depth_err2_idx = findfirst(x->x=="koi_depth_err2",csv_header)
  koi_duration_err2_idx = findfirst(x->x=="koi_duration_err2",csv_header)
  koi_disposition_idx = findfirst(x->x=="koi_disposition",csv_header)
  koi_pdisposition_idx = findfirst(x->x=="koi_pdisposition",csv_header)

  # TODO:  Let's move this subset selection code to a separate function, so we don't have so much duplicate code in this and previous function
  koi_subset = fill(false, length(csv_data[:,kepid_idx]))

  if haskey(sim_param, "koi_subset_csv")
    subset_df = readtable(convert(String,get(sim_param,"koi_subset_csv", "christiansen_kov.csv")), header=true, separator=' ')

    col_idx = findfirst(x->x==string(names(subset_df)[1]),csv_header)
    for n in 1:length(subset_df[:,1])
      subset_colnum = 1
      subset_entry = find(x->x==subset_df[n,1], csv_data[:,col_idx])
      # println("Initial cut: ", subset_entry)
      while (length(subset_entry) > 1) & (subset_colnum < length(names(subset_df)))
        subset_colnum += 1
        col_idx = findfirst(x->x==string(names(subset_df)[subset_colnum]),csv_header)
        subsubset = find(x->round(x*10.)==round(subset_df[n,subset_colnum]*10.), csv_data[subset_entry,col_idx])
	# println("Extra cut: ", subset_df[n,subset_colnum], " / ", csv_data[subset_entry,col_idx], " = ", subsubset)
	subset_entry = subset_entry[subsubset]
      end
      if subset_colnum > 1
 	col_idx = findfirst(x->x==string(names(subset_df)[1]),csv_header)
      end
      if length(subset_entry) > 1
	cand_sub = find(x->x == "CANDIDATE",csv_data[subset_entry,koi_pdisposition_idx])
	subset_entry = subset_entry[cand_sub]
	if length(subset_entry) > 1
          println("Multiple planets found in final cut: ", subset_df[n,1])
        end
      end
      if length(subset_entry) < 1
        println("No planets found in final cut: ", subset_df[n,:])
      end
      koi_subset[subset_entry] = true
    end
  else
    println("No KOI csv file given!")
    koi_subset = fill(true, length(csv_data[:,kepid_idx]))
  end

  csv_data = csv_data[koi_subset,:]
 
  output = KeplerObsCatalog([])

  tot_plan = count(x->x, koi_subset)
  for (j,kepid) in enumerate(df_star[:kepid])
     plids = find(x->x==kepid,csv_data[:,kepid_idx])
     num_pl = length(plids)
     tot_plan = tot_plan - num_pl
     if num_pl ==0 continue end
     target_obs = KeplerTargetObs(num_pl)
     target_obs.star = ExoplanetsSysSim.StarObs(df_star[j,:radius],df_star[j,:mass])
     for (plid,i) in enumerate(plids)
       target_obs.obs[plid] = ExoplanetsSysSim.TransitPlanetObs(csv_data[i,koi_period_idx],csv_data[i,koi_time0bk_idx],csv_data[i,koi_depth_idx]/1.0e6,csv_data[i,koi_duration_idx])
       target_obs.sigma[plid] = ExoplanetsSysSim.TransitPlanetObs((abs(csv_data[i,koi_period_err1_idx])+abs(csv_data[i,koi_period_err2_idx]))/2,(abs(csv_data[i,koi_time0bk_err1_idx])+abs(csv_data[i,koi_time0bk_err2_idx]))/2,(abs(csv_data[i,koi_depth_err1_idx]/1.0e6)+abs(csv_data[i,koi_depth_err2_idx]/1.0e6))/2,(abs(csv_data[i,koi_duration_err1_idx])+abs(csv_data[i,koi_duration_err2_idx]))/2)
       target_obs.prob_detect = ExoplanetsSysSim.SimulatedSystemDetectionProbs{OneObserver}( ones(num_pl), ones(num_pl,num_pl), ones(num_pl), fill(Array{Int64}(0), 1) )
       #target_obs.prob_detect = ExoplanetsSysSim.OneObserverSystemDetectionProbs(num_pl)
     end
     push!(output.target,target_obs)
  end
  println("Number of planet candidates with no matching star in table: ", tot_plan)
  return output
end


# Two functions below were just for debugging purposes
function calc_snr_list(cat::KeplerPhysicalCatalog, sim_param::SimParam)
  snrlist = Array{Float64}(0)
  for t in 1:length(cat.target)
    for p in 1:length(cat.target[t].sys[1].planet)
      snr = calc_snr_if_transit(cat.target[t],1,p,sim_param)
      if snr>0.0
        push!(snrlist,snr)
      end
    end
  end
  snrlist[find(x->x>7.1,snrlist)]
end

function calc_prob_detect_list(cat::KeplerPhysicalCatalog, sim_param::SimParam)
  pdetectlist = Array{Float64}(0)
  for t in 1:length(cat.target)
    for p in 1:length(cat.target[t].sys[1].planet)
      pdet = calc_prob_detect_if_transit(cat.target[t],1,p,sim_param)
      if pdet>0.0
        push!(pdetectlist,pdet)
      end
    end
  end
  idx = find(x->x>0.0,pdetectlist)
  pdetectlist[idx]
end

function test_catalog_constructors(sim_param::SimParam)
  cat_phys = generate_kepler_physical_catalog(sim_param)::KeplerPhysicalCatalog
  id = findfirst( x->num_planets(x)>=1 , cat_phys.target)   # fast forward to first target that has some planets
  @assert(length(id)>=1)
  semimajor_axis(cat_phys.target[id].sys[1],1)
  pdetlist = calc_prob_detect_list(cat_phys,sim_param)
  calc_target_obs_single_obs(cat_phys.target[id],sim_param)
  calc_target_obs_sky_ave(cat_phys.target[id],sim_param)
  @assert( length(cat_phys.target[id].sys[1].planet)  == num_planets(cat_phys.target[id]) )
  cat_obs = simulated_read_kepler_observations(sim_param)
  return (cat_phys, cat_obs)
end


