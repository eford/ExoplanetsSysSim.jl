## ExoplanetsSysSim/examples/hsu_etal_2018/calc_sum_stats_pr.jl
## 2017 Keir A. Ashby
#  Calculates summary statistics that can be used to generate period-radius distribution plots
#  Based primarily off of calc_summary_stats_obs_binned_rates from christiansen_func.jl

function calc_sum_stats_pr(cat_obs::KeplerObsCatalog, param::SimParam, trueobs_cat::Bool = false, cat_phys::KeplerPhysicalCatalog=KeplerPhysicalCatalog([]))
  ssd = Dict{ASCIIString,Any}()
  cache = Dict{ASCIIString,Any}()

  if !trueobs_cat
    num_targets = get_int(param,"num_targets_sim_pass_one")
    ssd["num targets"] = num_targets
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

    period_obs_list = zeros(num_tranets)
    weight_list = zeros(num_tranets)
    radius_obs_list = zeros(num_tranets)

  if trueobs_cat==true
    mult_true   = zeros(max_tranets_in_sys)

    n = 1    # tranet id
    for i in idx_tranets
      for j in 1:num_planets(cat_obs.target[i])
	if length(cat_obs.target[i].obs)==0 continue end	
        period_obs_list[n] = cat_obs.target[i].obs[j].period
        weight_list[n] = 1.0
        radius_obs_list[n] = sqrt(cat_obs.target[i].obs[j].depth)*cat_obs.target[i].star.radius
        num_planets_in_sys = length(cat_obs.target[i].obs)
        mult_true[num_planets_in_sys]+=1
        n = n+1
      end
    end
  
    ssd["period_list"] = period_obs_list
    ssd["radius_list"] = radius_obs_list
    ssd["multiplicity_true"] = mult_true

  else  #Here we'll calculate periods and radii for the physical and observed catalogs as well as their planet multiplicities.
    period_phys_list = Float64[] # zeros(length(idx_phys))      #Holds the periods for every planet in the physical catalog
    radius_phys_list = Float64[] # zeros(length(idx_phys))	#Holds the radii for every planet in the physical catalog
    mult_phys        = zeros(max_tranets_in_sys) 		#Holds the number of systems for each variation of planet multiplicity in the physical catalog
    mult_obs         = zeros(max_tranets_in_sys) 		#Holds the number of systems for each variation of planet multiplicity in the observed catalog     
    a = 1    # tranet id
    b = 1    # tranet id
    for i in 1:length(cat_phys.target) 				#get from physical catalog
       for s in 1:length(cat_phys.target[i].sys)
         if length(cat_phys.target[i].sys)==0 continue end
	 for p in 1:length(cat_phys.target[i].sys[s].planet)
	   if length(cat_phys.target[i].sys[s].planet)==0 continue end	
	   #period_phys_list[a] = cat_phys.target[i].sys[s].orbit[p].P
           #radius_phys_list[a] = cat_phys.target[i].sys[s].planet[p].radius
	   push!(period_phys_list,cat_phys.target[i].sys[s].orbit[p].P)
           push!(radius_phys_list,cat_phys.target[i].sys[s].planet[p].radius)
	   num_planets_in_sys = length(cat_phys.target[i].sys[s].planet)
           mult_phys[num_planets_in_sys]+=1
           a+=1
         end	
	end
    end
    for j in idx_tranets
      for h in 1:num_planets(cat_obs.target[j])
	if length(cat_obs.target[j].obs)==0 continue end	
        period_obs_list[b] = cat_obs.target[j].obs[h].period
	weight_list[b] = 1.0
        radius_obs_list[b] = sqrt(cat_obs.target[j].obs[h].depth*cat_obs.target[j].star.radius)
        num_planets_in_sys = length(cat_obs.target[j].obs)
        mult_obs[num_planets_in_sys]+=1
        b+=1
      end 
    end

    ssd["period_phys_list"]   = period_phys_list
    ssd["radius_phys_list"]   = radius_phys_list
    ssd["period_obs_list"]    = period_obs_list
    ssd["radius_obs_list"]    = radius_obs_list
    ssd["max_tranets_in_sys"] = max_tranets_in_sys
    ssd["multiplicity_phys"]  = mult_phys
    ssd["multiplicity_obs"]   = mult_obs
  end

  limitP::Array{Float64,1} = get_any(param, "p_lim_arr", Array{Float64,1})
  limitRp::Array{Float64,1} = get_any(param, "r_lim_arr", Array{Float64,1})

  np_bin = zeros((length(limitP)-1) * (length(limitRp)-1))
  np_bin_idx = 1
  for i in 1:(length(limitP)-1)
    P_match = find(x -> ((x > limitP[i]) && (x < limitP[i+1])), period_obs_list)
    for j in 1:(length(limitRp)-1)
      R_match = find(x -> ((x > limitRp[j]) && (x < limitRp[j+1])), radius_obs_list)
      
      bin_match = intersect(P_match, R_match)

      np_bin[np_bin_idx] = sum(weight_list[bin_match])
      np_bin_idx += 1
    end
  end

  #ssd["planets detected"] = sum(np_bin)
  ssd["planets table"] = np_bin

  return CatalogSummaryStatistics(ssd, cache)
end

