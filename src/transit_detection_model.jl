## ExoplanetsSysSim/src/transit_detection_model.jl
## (c) 2015 Eric B. Ford

# Several functions below based on https://github.com/christopherburke/KeplerPORTs/blob/master/KeplerPORTs_utils.py
# That follows the procedure outlined in Burke et al.(2015).
# However we don't currently interpolate the CDDP or mesthreshold to the relevant duration

function real_log_choose(m::Real, n::Real)
  lgamma(m+1)-lgamma(n+1)-lgamma(m-n+1.0)
end
  
function real_binom(k::Real, BigM::Real, f::Real)
    F1 = real_log_choose(BigM,k)
    F2 = k*log(f)
    F3 = (BigM-k)*log(1.0-f)
    x = exp(F1+F2+F3)
    return x
end
  
function kepler_window_function(num_transits_no_gaps::Real, duty_cycle::Real; min_transits::Real = 3.0)
  if num_transits_no_gaps < min_transits
     return 0.0
  else
     return max(1.0 - real_binom(min_transits,num_transits_no_gaps,duty_cycle), 0.0)
  end
end
    
# See Christiansen et al. (2015)  This assumes a linear limbdarkening coefficient of 0.6
function frac_depth_to_tps_depth(frac_depth::Real)
    const alp = 1.0874
    const bet = 1.0187
    const REALDEPTH2TPSSQUARE = 1.0  # WARNING: Waiting for this to be confirmed
    k = sqrt(frac_depth)
    tps_depth = min( (alp-bet*k) * frac_depth* REALDEPTH2TPSSQUARE, 1.0)   # NOTE: I added the max based on common sense
    return tps_depth
end

function detection_efficiency_theory(mes::Real; min_pdet_nonzero::Float64 = 0.0)
   const muoffset =  0.0
   const sig =  1.0
   const mesthresh = 7.1
   if mes > (9.0 - mesthresh - muoffset)
      return 1.0
   else
      pdet = 0.5 + 0.5*erf((mes - mesthresh - muoffset) / sqrt(2.0*sig*sig))
      pdet = pdet >= min_pdet_nonzero ? pdet : 0.0
      return pdet
   end
end

function detection_efficiency_fressin2013(mes::Real)
   const mesmin =  6.0
   const mesmax =  16.0
   if mes <= mesmin
      return 0.0
   elseif mes >= mesmax
      return 1.0
   else
     return (mes - mesmin) / (mesmax - mesmin)
  end
end


function detection_efficiency_christiansen2015(mes::Real; mes_threshold::Real = 7.1, min_pdet_nonzero::Float64 = 0.0)
   const a =  4.65  # from code for detection_efficiency(...) at https://github.com/christopherburke/KeplerPORTs/blob/master/KeplerPORTs_utils.py
   const b =  0.98
   #const a =  4.35 # from arxiv abstract. TODO DETAIL: Figure out which is better.  Informal testing showed it didn't matter
   #const b =  1.05
   usemes = max(0.0,mes - 7.1 - (mes_threshold - 7.1))
   pdet = cdf(Gamma(a,b), usemes)
   pdet = pdet >= min_pdet_nonzero ? pdet : 0.0
   return pdet
end

function detection_efficiency_dr25_simple(mes::Real; min_pdet_nonzero::Float64 = 0.0)
   const a =  30.87  # from pg 16 of https://exoplanetarchive.ipac.caltech.edu/docs/KSCI-19110-001.pdf
   const b =  0.271
   const c = 0.940
   pdet = c*cdf(Gamma(a,b), mes)
   pdet = pdet >= min_pdet_nonzero ? pdet : 0.0
   return pdet
end

detection_efficiency_model = detection_efficiency_christiansen2015  #  WARNING: Hardcoded choice of transit detection efficiency here for speed and so as to not have it hardcoded in multiple places

# Resume code original to SysSim

function calc_snr_if_transit(t::KeplerTarget, depth::Real, duration::Real, sim_param::SimParam; num_transit::Real = 1)
  cdpp = t.cdpp[1,1]                                                                # TODO: Make CDPP lookup based on season/quarter/month and timescale.  Is CDPP timescale IMPORTANT or DETAIL?
  depth_tps = frac_depth_to_tps_depth(depth)                                        # WARNING: Hardcoded this conversion
  snr = depth_tps*sqrt(num_transit*duration*LC_rate)/cdpp              # WARNING: Assumes measurement uncertainties are uncorrelated & CDPP based on LC
end

function calc_snr_if_transit(t::KeplerTarget, s::Integer, p::Integer, sim_param::SimParam)
  depth = calc_transit_depth(t,s,p)
  duration = calc_transit_duration(t,s,p)
  num_transit = calc_expected_num_transits(t,s,p,sim_param)
  calc_snr_if_transit(t,depth,duration,sim_param,num_transit=num_transit)
end

function calc_prob_detect_if_transit(t::KeplerTarget, snr::Real, sim_param::SimParam; num_transit::Real = 1)
  const min_transits = 3.0                                                    # WARNING: Hard coded 3 transit minimum
  const mes_threshold = 7.1                                                   # WARNING: Assuming 7.1 for all stars, durations
  const min_pdet_nonzero = 0.0                                                # TODO OPT: Figure out how to prevent a plethora of planets that are very unlikely to be detected due to using 0.0
  wf = kepler_window_function(num_transit, t.duty_cycle, min_transits=min_transits)   # TODO SCI DETAIL: Replace statistical model with checking actual transit times for long period planets
  return wf*detection_efficiency_model(snr, min_pdet_nonzero=min_pdet_nonzero)	
end

function calc_prob_detect_if_transit(t::KeplerTarget, depth::Real, duration::Real, sim_param::SimParam; num_transit::Real = 1)
  snr = calc_snr_if_transit(t,depth,duration,sim_param, num_transit=num_transit)
  return calc_prob_detect_if_transit(t, snr, sim_param, num_transit=num_transit)
end

function calc_prob_detect_if_transit(t::KeplerTarget, s::Integer, p::Integer, sim_param::SimParam)
  depth = calc_transit_depth(t,s,p)
  duration = calc_transit_duration(t,s,p)
  ntr = calc_expected_num_transits(t,s,p,sim_param)
  calc_prob_detect_if_transit(t,depth,duration, sim_param, num_transit=ntr)
end

# Compute probability of detection if we average over impact parameters b~U[0,1)
function calc_ave_prob_detect_if_transit(t::KeplerTarget, snr_central::Real, sim_param::SimParam; num_transit::Real = 1)
  const min_transits = 3.0                                                    # WARNING: Hard coded 3 transit minimum
  const mes_threshold = 7.1                                                   # WARNING: Assuming 7.1 for all stars, durations
  const min_pdet_nonzero = 0.0                                                # TODO OPT: Figure out how to prevent a plethora of planets that are very unlikely to be detected due to using 0.0
  wf = kepler_window_function(num_transit, t.duty_cycle, min_transits=min_transits)   
  num_impact_param = 5
  ave_detection_efficiency = 0.5*detection_efficiency_model(snr_central, min_pdet_nonzero=min_pdet_nonzero)	  
  for i in 1:(num_impact_param-1)
    b = i/num_impact_param
    snr = snr_central*sqrt(sqrt((1-b)*(1+b)))
    ave_detection_efficiency += detection_efficiency_model(snr, min_pdet_nonzero=min_pdet_nonzero)
  end
  ave_detection_efficiency /= num_impact_param
  return wf*ave_detection_efficiency
end

function calc_ave_prob_detect_if_transit(t::KeplerTarget, depth::Real, duration_central::Real, sim_param::SimParam; num_transit::Real = 1)
  snr_central = calc_snr_if_transit(t,depth,duration_central,sim_param, num_transit=num_transit)
  return calc_ave_prob_detect_if_transit(t, snr_central, sim_param, num_transit=num_transit)
end

function calc_ave_prob_detect_if_transit(t::KeplerTarget, s::Integer, p::Integer, sim_param::SimParam)
  depth = calc_transit_depth(t,s,p)
  duration_central = calc_transit_duration_central(t,s,p)
  ntr = calc_expected_num_transits(t,s,p,sim_param)
  calc_ave_prob_detect_if_transit(t,depth,duration_central, sim_param, num_transit=ntr)
end
