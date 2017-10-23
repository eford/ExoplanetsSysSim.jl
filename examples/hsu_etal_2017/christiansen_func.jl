using ExoplanetsSysSim
using StatsFuns
using JLD
using DataFrames
#import ExoplanetsSysSim.StellarTable.df
#import ExoplanetsSysSim.StellarTable.usable
#import Compat: UTF8String, ASCIIString

## simulation_parameters
function setup_sim_param_christiansen(args::Vector{String} = Array{String}(0) )   # allow this to take a list of parameter (e.g., from command line)
  sim_param = setup_sim_param_demo()
  
  set_inactive(sim_param,["log_eta_pl","power_law_P","power_law_r"])
  add_param_fixed(sim_param,"star_table_setup",setup_star_table_christiansen)
  add_param_fixed(sim_param,"stellar_catalog","q1_q17_christiansen.jld")
  add_param_fixed(sim_param,"generate_num_planets",generate_num_planets_christiansen)
  add_param_fixed(sim_param,"generate_period_and_sizes", generate_period_and_sizes_christiansen)
  add_param_fixed(sim_param,"num_targets_sim_pass_one",150000)  # For faster simulated catalogs

  p_lim_arr_num = [0.5, 1.25, 2.5, 5., 10., 20., 40., 80., 160., 320.]
  r_lim_arr_num = [1., 1.25, 1.5, 1.75, 2.]
  rate_tab_init = [0.027 0.149 0.354 1.37 2.56 1.96 2.85 4.77 0.0001 ;
                   0.02 0.063 0.432 0.9 1.72 1.46 2.4 5.42 3.63 ;
                   0.0074 0.032 0.389 0.545 1.32 1.65 1.96 1.69 1.17 ;
                   0.0001 0.017 0.156 0.542 0.777 1.35 0.965 1.22 2.49]*0.01

  add_param_fixed(sim_param, "p_lim_arr", p_lim_arr_num)
  add_param_fixed(sim_param, "r_lim_arr", r_lim_arr_num*ExoplanetsSysSim.earth_radius)
  add_param_active(sim_param,"obs_par", rate_tab_init)

  return sim_param
end


## planetary_system
function generate_num_planets_christiansen(s::Star, sim_param::SimParam)
  const max_tranets_in_sys::Int64 = get_int(sim_param,"max_tranets_in_sys")
  rate_tab::Array{Float64,2} = get_any(sim_param, "obs_par", Array{Float64,2})
  lambda = sum_kbn(rate_tab)
  #println("# lambda= ", lambda)
  ExoplanetsSysSim.generate_num_planets_poisson(lambda,max_tranets_in_sys)
end

function generate_period_and_sizes_christiansen(s::Star, sim_param::SimParam; num_pl::Integer = 1)
  rate_tab::Array{Float64,2} = get_any(sim_param, "obs_par", Array{Float64,2})
  
  limitP::Array{Float64,1} = get_any(sim_param, "p_lim_arr", Array{Float64,1})
  limitRp::Array{Float64,1} = get_any(sim_param, "r_lim_arr", Array{Float64,1})

  @assert ((length(limitP)-1) == size(rate_tab, 2))
  @assert ((length(limitRp)-1) == size(rate_tab, 1))
  Plist = []
  Rplist = []
  rate_tab_1d = reshape(rate_tab,length(rate_tab))
  #logmaxcuml = logsumexp(rate_tab_1d)
  #cuml = cumsum_kbn(exp(rate_tab_1d-logmaxcuml))
  maxcuml = sum(rate_tab_1d)
  cuml = cumsum_kbn(rate_tab_1d/maxcuml)

  for n in 1:num_pl
    rollp = rand()
    idx = findfirst(x -> x > rollp, cuml)
    i_idx = (idx-1)%size(rate_tab,1)+1
    j_idx = floor(Int64,(idx-1)//size(rate_tab,1))+1
    ### TODO: Keep uniform log sampling here?
    Rp = exp(rand()*(log(limitRp[i_idx+1])-log(limitRp[i_idx]))+log(limitRp[i_idx]))
    P = exp(rand()*(log(limitP[j_idx+1])-log(limitP[j_idx]))+log(limitP[j_idx]))
    push!(Plist, P)
    push!(Rplist, Rp)
  end
  return Plist, Rplist
end


## stellar_table
function setup_christiansen(sim_param::SimParam; force_reread::Bool = false)
  #global df
  df = ExoplanetsSysSim.StellarTable.df
  if haskey(sim_param,"read_stellar_catalog") && !force_reread
     return df
     #return data
  end
  stellar_catalog_filename = convert(String,joinpath(Pkg.dir("ExoplanetsSysSim"), "data", convert(String,get(sim_param,"stellar_catalog","q1_q17_dr25_stellar.csv")) ) )
  df = setup_christiansen(stellar_catalog_filename)
  add_param_fixed(sim_param,"read_stellar_catalog",true)
  set_star_table(df)
  return df  
end

function setup_christiansen(filename::String; force_reread::Bool = false)
  #global df, usable
  df = ExoplanetsSysSim.StellarTable.df
  usable = ExoplanetsSysSim.StellarTable.usable
  if ismatch(r".jld$",filename)
  try 
    data = load(filename)
    df::DataFrame = data["stellar_catalog"]
    usable::Array{Int64,1} = data["stellar_catalog_usable"]
    set_star_table(df, usable)
  catch
    error(string("# Failed to read stellar catalog >",filename,"< in jld format."))
  end
  else
  try 
    df = readtable(filename)
  catch
    error(string("# Failed to read stellar catalog >",filename,"< in ascii format."))
  end

  has_mass = ! (isna(df[:mass]) | isna(df[:mass_err1]) | isna(df[:mass_err2]))
  has_radius = ! (isna(df[:radius]) | isna(df[:radius_err1]) | isna(df[:radius_err2]))
  has_dens = ! (isna(df[:dens]) | isna(df[:dens_err1]) | isna(df[:dens_err2]))
  has_rest = ! (isna(df[:rrmscdpp04p5]) | isna(df[:dataspan]) | isna(df[:dutycycle]))
  in_Q1Q12 = []
  for x in df[:st_quarters]
    subx = string(x)
    subx = ("0"^(17-length(subx)))*subx
    indQ = search(subx, '1')
    if ((indQ < 1) | (indQ > 12))
      push!(in_Q1Q12, false)
    else
      push!(in_Q1Q12, true)
    end
  end
  is_FGK = []
  for x in 1:length(df[:teff])
    if ((df[x,:teff] > 4000.0) & (df[x,:teff] < 7000.0) & (df[x,:logg] > 4.0))
      push!(is_FGK, true)
    else
      push!(is_FGK, false)
    end
  end
  is_usable = has_radius & is_FGK & has_mass & has_rest #& in_Q1Q12 # & has_dens
  if contains(filename,"q1_q12_christiansen.jld")
    is_usable = is_usable & in_Q1Q12
  end
  # See options at: http://exoplanetarchive.ipac.caltech.edu/docs/API_keplerstellar_columns.html
  # TODO SCI DETAIL or IMPORTANT?: Read in all CDPP's, so can interpolate?
  symbols_to_keep = [ :kepid, :mass, :mass_err1, :mass_err2, :radius, :radius_err1, :radius_err2, :dens, :dens_err1, :dens_err2, :rrmscdpp04p5, :dataspan, :dutycycle ]
  delete!(df, [~(x in symbols_to_keep) for x in names(df)])    # delete columns that we won't be using anyway
  usable = find(is_usable)
  df = df[usable, symbols_to_keep]
  set_star_table(df, usable)
  end
  return df
  #global data = convert(Array{Float64,2}, df) # df[usable, symbols_to_keep] )
  #global colid = Dict(zip(names(df),[1:length(names(df))]))
  #return data
end

setup_star_table_christiansen(sim_param::SimParam; force_reread::Bool = false) = setup_christiansen(sim_param, force_reread=force_reread)
setup_star_table_christiansen(filename::String; force_reread::Bool = false) = setup_christiansen(filename, force_reread=force_reread)


## summary_statistics
function calc_summary_stats_sim_pass_one_binned_rates(cat_obs::KeplerObsCatalog, cat_phys::KeplerPhysicalCatalog, param::SimParam )      # Version for simulated data, since includes cat_phys
  ssd = Dict{String,Any}()
  cache = Dict{String,Any}()

  max_tranets_in_sys = get_int(param,"max_tranets_in_sys")    # Demo that simulation parameters can specify how to evalute models, too
  @assert max_tranets_in_sys >= 1
  idx_tranets = find(x::KeplerTargetObs-> length(x.obs) > 0, cat_obs.target)::Array{Int64,1}             # Find indices of systems with at least 1 tranet = potentially detectable transiting planet

  # Count total number of tranets and compile indices for N-tranet systems
  num_tranets = 0
  idx_n_tranets = Vector{Int64}[ Int64[] for m = 1:max_tranets_in_sys]
  for n in 1:max_tranets_in_sys-1
    idx_n_tranets[n] = find(x::KeplerTargetObs-> length(x.obs) == n, cat_obs.target[idx_tranets] )
    num_tranets += n*length(idx_n_tranets[n])
  end
  idx_n_tranets[max_tranets_in_sys] = find(x::KeplerTargetObs-> length(x.obs) >= max_tranets_in_sys, cat_obs.target[idx_tranets] )

  num_tranets += max_tranets_in_sys*length(idx_n_tranets[max_tranets_in_sys])  # WARNING: this means we need to ignore planets w/ indices > max_tranets_in_sys
  num_tranets  = convert(Int64,num_tranets)            # TODO OPT: Figure out why isn't this already an Int.  I may be doing something that prevents some optimizations

  cache["num_tranets"] = num_tranets                                   
  cache["idx_tranets"] = idx_tranets                                   # We can save lists of indices to summary stats for pass 2, even though we won't use these for computing a distance or probability
  #cache["idx_n_tranets"] = idx_n_tranets

  expected_num_detect = 0.0
  expected_num_sys_n_tranets = zeros(max_tranets_in_sys)
  period_list = zeros(num_tranets)
  weight_list = zeros(num_tranets)
  radius_list = zeros(num_tranets)

  n = 1    # tranet id
  for i in idx_tranets
    for j in 1:num_planets(cat_obs.target[i])
      p_tr_and_det = ExoplanetsSysSim.prob_detect(cat_obs.target[i].prob_detect,j)
      expected_num_detect += p_tr_and_det
      (s,p) = cat_obs.target[i].phys_id[j]
      
      period_list[n] = cat_phys.target[i].sys[s].orbit[p].P
      weight_list[n] = p_tr_and_det
      radius_list[n] = cat_phys.target[i].sys[s].planet[p].radius
      n = n+1
    end
    for k in 1:max_tranets_in_sys
      expected_num_sys_n_tranets[k] += ExoplanetsSysSim.prob_detect_n_planets(cat_obs.target[i].prob_detect,k)
    end
  end
  ssd["expected planets detected"] = expected_num_detect
  ssd["num_sys_tranets"] = expected_num_sys_n_tranets
  ssd["num targets"] = get_int(param,"num_targets_sim_pass_one")
  #println("expected planets = ",expected_num_detect,", num_sys_tranets = ",expected_num_sys_n_tranets,", num targets = ",ssd["num targets"])

  limitP::Array{Float64,1} = get_any(param, "p_lim_arr", Array{Float64,1})
  limitRp::Array{Float64,1} = get_any(param, "r_lim_arr", Array{Float64,1})

  np_bin = zeros((length(limitP)-1) * (length(limitRp)-1))
  np_bin_idx = 1
  for i in 1:(length(limitP)-1)
    P_match = find(x -> ((x > limitP[i]) && (x < limitP[i+1])), period_list)
    for j in 1:(length(limitRp)-1)
      R_match = find(x -> ((x > limitRp[j]) && (x < limitRp[j+1])), radius_list)
      
      bin_match = intersect(P_match, R_match)

      np_bin[np_bin_idx] = sum(weight_list[bin_match])
      np_bin_idx += 1
    end
  end

  #ssd["expected planets detected"] = sum(np_bin)
  ssd["expected planets table"] = np_bin

  return CatalogSummaryStatistics(ssd, cache)
end

function calc_summary_stats_obs_binned_rates(cat_obs::KeplerObsCatalog, param::SimParam; trueobs_cat::Bool = false)
  ssd = Dict{String,Any}()
  cache = Dict{String,Any}()

  if !trueobs_cat
    ssd["num targets"] = get_int(param,"num_targets_sim_pass_one")
  else
    ssd["num targets"] = get_int(param,"num_kepler_targets")
  end

  max_tranets_in_sys = get_int(param,"max_tranets_in_sys")    # Demo that simulation parameters can specify how to evalute models, too
  @assert max_tranets_in_sys >= 1
  idx_tranets = find(x::KeplerTargetObs-> length(x.obs) > 0, cat_obs.target)::Array{Int64,1}             # Find indices of systems with at least 1 tranet = potentially detectable transiting planet

  # Count total number of tranets and compile indices for N-tranet systems
  num_tranets = 0
  idx_n_tranets = Vector{Int64}[ Int64[] for m = 1:max_tranets_in_sys]
  for n in 1:max_tranets_in_sys-1
    idx_n_tranets[n] = find(x::KeplerTargetObs-> length(x.obs) == n, cat_obs.target[idx_tranets] )
    num_tranets += n*length(idx_n_tranets[n])
  end
  idx_n_tranets[max_tranets_in_sys] = find(x::KeplerTargetObs-> length(x.obs) >= max_tranets_in_sys, cat_obs.target[idx_tranets] )

  num_tranets += max_tranets_in_sys*length(idx_n_tranets[max_tranets_in_sys])  # WARNING: this means we need to ignore planets w/ indices > max_tranets_in_sys
  if ( length( find(x::KeplerTargetObs-> length(x.obs) > max_tranets_in_sys, cat_obs.target[idx_tranets] ) ) > 0)   # Make sure max_tranets_in_sys is at least big enough for observed systems
    warn("Observational data has more transiting planets in one systems than max_tranets_in_sys allows.")
  end
  num_tranets  = convert(Int64,num_tranets)            # TODO OPT: Figure out why isn't this already an Int.  I may be doing something that prevents some optimizations

  num_sys_tranets = zeros(max_tranets_in_sys)                           # Since observed data, don't need to calculate probabilities.
  for n in 1:max_tranets_in_sys                                         # Make histogram of N-tranet systems
    num_sys_tranets[n] = length(idx_n_tranets[n])
  end
  ssd["num_sys_tranets"] = num_sys_tranets
  ssd["planets detected"] = num_tranets 

  period_list = zeros(num_tranets)
  weight_list = zeros(num_tranets)
  radius_list = zeros(num_tranets)

  n = 1    # tranet id
  for i in idx_tranets
    for j in 1:num_planets(cat_obs.target[i])
      period_list[n] = cat_obs.target[i].obs[j].period
      weight_list[n] = 1.0
      radius_list[n] = sqrt(cat_obs.target[i].obs[j].depth)*cat_obs.target[i].star.radius
      n = n+1
    end
  end

  limitP::Array{Float64,1} = get_any(param, "p_lim_arr", Array{Float64,1})
  limitRp::Array{Float64,1} = get_any(param, "r_lim_arr", Array{Float64,1})

  np_bin = zeros((length(limitP)-1) * (length(limitRp)-1))
  np_bin_idx = 1
  for i in 1:(length(limitP)-1)
    P_match = find(x -> ((x > limitP[i]) && (x < limitP[i+1])), period_list)
    for j in 1:(length(limitRp)-1)
      R_match = find(x -> ((x > limitRp[j]) && (x < limitRp[j+1])), radius_list)
      
      bin_match = intersect(P_match, R_match)

      np_bin[np_bin_idx] = sum(weight_list[bin_match])
      np_bin_idx += 1
    end
  end

  #ssd["planets detected"] = sum(np_bin)
  ssd["planets table"] = np_bin

  return CatalogSummaryStatistics(ssd, cache)
end


## abc_distance
function calc_distance_vector_binned(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, pass::Int64, sim_param::SimParam ; verbose::Bool = false)
  d = Array{Float64}(0)
  if pass == 1
    if verbose
      println("# Summary 1, pass 1: ",summary1)
      println("# Summary 2, pass 1: ",summary2)
    end
    d = zeros(3)
    # Since observed and simulated catalogs can have different summary statistics for the number of planets, prefer detections if avaliable (e.g., after pass2), otherwise use expected (e.g., from pass 1)
    #np1 = haskey(summary1.stat,"planets detected") ? summary1.stat["planets detected"] : summary1.stat["expected planets detected"]
    #np2 = haskey(summary2.stat,"planets detected") ? summary2.stat["planets detected"] : summary2.stat["expected planets detected"]
    #d[1] = dist_L1_abs(np1/summary1.stat["num targets"],np2/summary2.stat["num targets"])    #  Normalize so different statistics weighted appropriately and not dominated by this one
    #println("np1 = ",np1,", np2 = ",np2)
    #println("np1 (normalized) = ",np1/summary1.stat["num targets"],", np2 (normalized) = ",np2/summary2.stat["num targets"],", d[1] = ",d[1])

    np1 = haskey(summary1.stat,"planets table") ? summary1.stat["planets table"] : summary1.stat["expected planets table"]
    np2 = haskey(summary2.stat,"planets table") ? summary2.stat["planets table"] : summary2.stat["expected planets table"]

    np_bin = zeros(length(np1))
    for n in 1:length(np1)
        #np_bin[n] = dist_L1_abs(np1[n]/summary1.stat["num targets"], np2[n]/summary2.stat["num targets"])
        np_bin[n] = dist_L2_abs(np1[n]/summary1.stat["num targets"], np2[n]/summary2.stat["num targets"])
      #println("True # [Bin ", n,"] = ",np1[n],", Expected # [Bin ", n,"] = ",np2[n])
    end
      #d[1] = maximum(np_bin)
      d[1] = sum(np_bin)
    else
    println("# calc_distance_vector_demo doesn't know what to do for pass= ", pass)
  end
  return d
end


## eval_model
function test_christiansen()
  global sim_param_closure = setup_sim_param_christiansen()
  cat_phys = generate_kepler_physical_catalog(sim_param_closure)
  cat_obs = observe_kepler_targets_single_obs(cat_phys,sim_param_closure)
  global summary_stat_ref_closure = calc_summary_stats_obs_demo(cat_obs,sim_param_closure)
  global cat_phys_try_closure  = generate_christiansen_catalog(sim_param_closure)
  global cat_obs_try_closure  = observe_kepler_targets_sky_avg(cat_phys_try_closure,sim_param_closure)
  global summary_stat_try_closure  = calc_summary_stats_sim_pass_one_demo(cat_obs_try_closure,cat_phys_try_closure,sim_param_closure)
  summary_stat_try_closure   = calc_summary_stats_sim_pass_two_demo(cat_obs_try_closure,cat_phys_try_closure,summary_stat_try_closure,sim_param_closure)
  param_guess = make_vector_of_sim_param(sim_xparam_closure)
  evaluate_model_scalar_ret( param_guess)
end


## inverse_detection
function inv_det_prob(param::SimParam)
  cat_obs = setup_actual_planet_candidate_catalog(setup_star_table_christiansen(param), param)
  inv_det_prob(cat_obs, param)
end

function inv_det_prob(cat_obs::KeplerObsCatalog, param::SimParam)
  num_targ = num_usable_in_star_table()
  idx_tranets = find(x::KeplerTargetObs-> length(x.obs) > 0, cat_obs.target)::Array{Int64,1} 

  limitP::Array{Float64,1} = get_any(param, "p_lim_arr", Array{Float64,1})
  limitRp::Array{Float64,1} = get_any(param, "r_lim_arr", Array{Float64,1})

  np_bin = zeros((length(limitP)-1) * (length(limitRp)-1))
  w2_bin = zeros((length(limitP)-1) * (length(limitRp)-1))
  cnt_bin = zeros((length(limitP)-1) * (length(limitRp)-1))
  pl_idx = 1

  for i in idx_tranets
    for j in 1:num_planets(cat_obs.target[i])
      pper = cat_obs.target[i].obs[j].period
      prad = sqrt(cat_obs.target[i].obs[j].depth)*cat_obs.target[i].star.radius

      pbin = findfirst(x -> ((pper > limitP[x]) && (pper < limitP[x+1])), collect(1:(length(limitP)-1)))
      rbin = findfirst(x -> ((prad > limitRp[x]) && (prad < limitRp[x+1])), collect(1:(length(limitRp)-1)))
      if (pbin > 0 && rbin > 0)
        cnt_bin[(pbin-1)*(length(limitRp)-1) + rbin] += 1
        pgeo = ExoplanetsSysSim.calc_transit_prob_single_planet_approx(pper, cat_obs.target[i].star.radius, cat_obs.target[i].star.mass)
	pdet = 0.0
	for star_id in 1:num_targ
	  star = SingleStar(star_table(star_id,:radius),star_table(star_id,:mass),1.0, star_id)
	  cdpp = 1.0e-6 * star_table(star_id, :rrmscdpp04p5) * sqrt(4.5/24.0 / ExoplanetsSysSim.LC_duration )
	  contam = 0.0
	  data_span = star_table(star_id, :dataspan)
	  duty_cycle = star_table(star_id, :dutycycle)
	  pl_arr = Array{Planet}( 1)
	  orbit_arr = Array{Orbit}( 1)
          incl = acos(rand()*star.radius*ExoplanetsSysSim.rsol_in_au/ExoplanetsSysSim.semimajor_axis(pper, star.mass))
	  orbit_arr[1] = Orbit(pper, 0., incl, 0., 0., rand()*2.*pi)
	  pl_arr[1] = Planet(prad, 1.0e-6)
	  kep_targ = KeplerTarget([PlanetarySystem(star, pl_arr, orbit_arr)], fill(cdpp,ExoplanetsSysSim.num_cdpp_timescales,ExoplanetsSysSim.num_quarters),contam,data_span,duty_cycle)

	  duration_central = ExoplanetsSysSim.calc_transit_duration(kep_targ,1,1) 
	  if duration_central <= 0.
	    continue
	  end
	  ntr = ExoplanetsSysSim.calc_expected_num_transits(kep_targ, 1, 1, param)
	  depth = ExoplanetsSysSim.calc_transit_depth(kep_targ,1,1)
          snr_central = ExoplanetsSysSim.calc_snr_if_transit(kep_targ, depth, duration_central, param, num_transit=ntr)
	  pdet += ExoplanetsSysSim.calc_prob_detect_if_transit(kep_targ, snr_central, param, num_transit=ntr)
	end
        np_bin[(pbin-1)*(length(limitRp)-1) + rbin] += 1.0/pgeo/(pdet/num_targ)
	w2_bin[(pbin-1)*(length(limitRp)-1) + rbin] += (1.0/pgeo/(pdet/num_targ))^2
	println("Planet ",pl_idx," => Bin ", (pbin-1)*(length(limitRp)-1) + rbin, ", C = ", 1.0/pgeo/(pdet/num_targ))
	pl_idx += 1
      end
    end
  end

  println()
  ess_bin = stellar_ess(param)
  println()

  println("Christiansen Rates:")
  for i in 1:(length(limitP)-1)
    for j in 1:(length(limitRp)-1)
      rate_f = np_bin[(i-1)*(length(limitRp)-1) + j]/num_targ*100.
      if cnt_bin[(i-1)*(length(limitRp)-1) + j] > 0.
        println(rate_f, 
                " +/- ", rate_f/sqrt(cnt_bin[(i-1)*(length(limitRp)-1) + j]), " %")
      else
        println(rate_f, 
                " +/- N/A %")
      end
    end
  end

  println()
  println("New Rates:")
  for i in 1:(length(limitP)-1)
    for j in 1:(length(limitRp)-1)
      rate_f = (1.0+cnt_bin[(i-1)*(length(limitRp)-1) + j])/(1.0+ess_bin[(i-1)*(length(limitRp)-1) + j])*100.
      up_quant = quantile(Gamma(1.0+cnt_bin[(i-1)*(length(limitRp)-1) + j], 1.0/(1.0+ess_bin[(i-1)*(length(limitRp)-1) + j])), 0.8413)*100.
      low_quant = quantile(Gamma(1.0+cnt_bin[(i-1)*(length(limitRp)-1) + j], 1.0/(1.0+ess_bin[(i-1)*(length(limitRp)-1) + j])), 0.1587)*100.
#      println(rate_f, 
#                " +/- ", rate_f/sqrt(1.+cnt_bin[(i-1)*(length(limitRp)-1) + j]), " %")
       println(rate_f, 
               " + ", up_quant - rate_f,
               " - ", rate_f - low_quant, " %")
    end
  end

end

## stellar catalog ess
function stellar_ess(param::SimParam)
  num_realiz = 100
  setup_star_table_christiansen(param)
  num_targ = num_usable_in_star_table()

  limitP::Array{Float64,1} = get_any(param, "p_lim_arr", Array{Float64,1})
  limitRp::Array{Float64,1} = get_any(param, "r_lim_arr", Array{Float64,1})

  np_bin = zeros((length(limitP)-1) * (length(limitRp)-1))
  w2_bin = zeros((length(limitP)-1) * (length(limitRp)-1))

  for star_id in 1:num_targ
    star = SingleStar(star_table(star_id,:radius),star_table(star_id,:mass),1.0, star_id)
    cdpp = 1.0e-6 * star_table(star_id, :rrmscdpp04p5) * sqrt(4.5/24.0 / ExoplanetsSysSim.LC_duration )
    contam = 0.0
    data_span = star_table(star_id, :dataspan)
    duty_cycle = star_table(star_id, :dutycycle)
  
    for i_idx in 1:(length(limitP)-1)
      for j_idx in 1:(length(limitRp)-1)
        temp_bin = 0.0
        for n_test in 1:num_realiz
          pper = exp(rand()*(log(limitP[i_idx+1])-log(limitP[i_idx]))+log(limitP[i_idx]))
	  prad = exp(rand()*(log(limitRp[j_idx+1])-log(limitRp[j_idx]))+log(limitRp[j_idx]))
    
          pgeo = ExoplanetsSysSim.calc_transit_prob_single_planet_approx(pper, star.radius, star.mass)
	  pdet = 0.0
	
	  pl_arr = Array{Planet}(1)
	  orbit_arr = Array{Orbit}(1)
          incl = acos(rand()*star.radius*ExoplanetsSysSim.rsol_in_au/ExoplanetsSysSim.semimajor_axis(pper, star.mass))
	  orbit_arr[1] = Orbit(pper, 0., incl, 0., 0., rand()*2.*pi)
	  pl_arr[1] = Planet(prad, 1.0e-6)
	  kep_targ = KeplerTarget([PlanetarySystem(star, pl_arr, orbit_arr)], fill(cdpp,ExoplanetsSysSim.num_cdpp_timescales,ExoplanetsSysSim.num_quarters),contam,data_span,duty_cycle)

	  duration_central = ExoplanetsSysSim.calc_transit_duration(kep_targ,1,1) 
	  if duration_central <= 0.
	    continue
	  end
	  ntr = ExoplanetsSysSim.calc_expected_num_transits(kep_targ, 1, 1, param)
	  depth = ExoplanetsSysSim.calc_transit_depth(kep_targ,1,1)
          snr_central = ExoplanetsSysSim.calc_snr_if_transit(kep_targ, depth, duration_central, param, num_transit=ntr)
	  #println(semimajor_axis(kep_targ.sys[1], 1))
	  #println(pper, " ", duration, " ", depth, " ", snr, " ", ntr)
	  pdet = ExoplanetsSysSim.calc_prob_detect_if_transit(kep_targ, snr_central, param, num_transit=ntr)
	
	  temp_bin += (pgeo*pdet)
          #np_bin[(i_idx-1)*(length(limitRp)-1) + j_idx] += (pgeo*pdet)
	  #w2_bin[(i_idx-1)*(length(limitRp)-1) + j_idx] += (pgeo*pdet)^2
        end
        np_bin[(i_idx-1)*(length(limitRp)-1) + j_idx] += temp_bin/num_realiz
      end
    end
    if (rem(star_id, 1000) == 0.)
      println(string("Star #", star_id, " finished"))
    end
  end
  
  #ess_bin = zeros((length(limitP)-1) * (length(limitRp)-1))

  for i in 1:(length(limitP)-1)
    for j in 1:(length(limitRp)-1)
  #    ess_bin[(i-1)*(length(limitRp)-1) + j] = (np_bin[(i-1)*(length(limitRp)-1) + j]^2)/w2_bin[(i-1)*(length(limitRp)-1) + j]
      println("Period limits: ", limitP[i:i+1], " / Radius limits: ", limitRp[j:j+1]/ExoplanetsSysSim.earth_radius, " / ESS = ", np_bin[(i-1)*(length(limitRp)-1) + j])
    end
  end
  return np_bin
end
