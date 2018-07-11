## ExoplanetsSysSim/src/target.jl
## (c) 2015 Eric B. Ford

#using Distributions

immutable KeplerTarget
  #sys::PlanetarySystem   # Make array for planetary systems aroud multiple stars in one target?
  sys::Vector{PlanetarySystemAbstract}
  cdpp::Array{Float64,2} # fractional, not ppm; 2D to allow for multiple time scales, months, quarters or seasons/spacecraft rotation angles
                           # QUERY: Do we want a separate CDPP for SC?  Or will CDPP's in different months be for LC/SC depending on this variable?
                            # QUERY: Should this be moved to KeplerTargetObs?
                            # QUERY: Should we not add this to target and use the star id to lookup CDPP from the stellar table?
  contam::Float64     # QUERY: Do we want/need this, since we're able to generate multiple stars in a single target?
  data_span::Float64
  duty_cycle::Float64
  window_function_id::Int64 # Points to the id of the window function for this target
  #channel::Int64         # E.g., if we cared which Kepler channel the target fell on
  #has_sc::Vector{Bool}   # TODO OPT: Make Immutable Vector or BitArray for speed?  QUERY: Should this go in KeplerTargetObs?
  #                       # QUERY: Do we want a separate CDPP for SC?  Or will CDPP's in different months be for LC/SC depending on this variable?
  #ra::Float64           # E.g., if we cared about position on sky     QUERY:  Should we replace with galactic longitude and latitute?
  #dec::Floa64           #
end
num_planets(t::KeplerTarget) = sum( num_planets, t.sys)
flux(t::KeplerTarget) = sum(flux,t.sys)+t.contam



function draw_asymmetric_normal(mu::Real, sig_plus::Real, sig_minus::Real)
  stdn = randn()
  mu + ( (stdn>=zero(stdn)) ? sig_plus*stdn : sig_minus*stdn )
end

function generate_kepler_target_from_table(sim_param::SimParam)  
  #const  generate_star = get_function(sim_param,"generate_star")
  const  generate_planetary_system = get_function(sim_param,"generate_planetary_system")
  const use_star_table_sigmas = false

  max_star_id = StellarTable.num_usable_in_star_table()

  star_table(id::Integer,sym::Symbol) = StellarTable.star_table(id,sym)

  @assert(1<=max_star_id)
  star_id = rand(1:max_star_id)
  mass = 0.0
  dens = 0.0
  radius = 0.0
  if use_star_table_sigmas
    while !(0.0<mass<100.0)
       mass = draw_asymmetric_normal( star_table(star_id,:mass), star_table(star_id,:mass_err1), star_table(star_id,:mass_err2) )
    end
    while !(0.0<dens<1000.0)
      dens = draw_asymmetric_normal( star_table(star_id,:dens), star_table(star_id,:dens_err1), star_table(star_id,:dens_err2) )
    end
    while !(0.0<radius<2.0)
      radius = draw_asymmetric_normal( star_table(star_id,:radius), star_table(star_id,:radius_err1), star_table(star_id,:radius_err2) )
    end
  else
    mass   = star_table(star_id,:mass)
    dens   = star_table(star_id,:dens)
    radius = star_table(star_id,:radius)
  end
  star = SingleStar(radius,mass,1.0, star_id)        # TODO SCI: Allow for blends, binaries, etc.
    #cdpp = 1.0e-6 * star_table(star_id, :rrmscdpp04p5) * sqrt(4.5/24.0 / LC_duration )  # TODO SCI: Allow for multiple timescales
  cdpp_arr = (1.0e-6*sqrt(1./24.0/LC_duration)) .* [star_table(star_id, :rrmscdpp01p5)*sqrt(1.5), star_table(star_id, :rrmscdpp02p0)*sqrt(2.), star_table(star_id,:rrmscdpp02p5)*sqrt(2.5), star_table(star_id,:rrmscdpp03p0)*sqrt(3.), star_table(star_id,:rrmscdpp03p5)*sqrt(3.5), star_table(star_id,:rrmscdpp04p5)*sqrt(4.5), star_table(star_id,:rrmscdpp05p0)*sqrt(5.), star_table(star_id,:rrmscdpp06p0)*sqrt(6.), star_table(star_id,:rrmscdpp07p5)*sqrt(7.5), star_table(star_id,:rrmscdpp09p0)*sqrt(9.), star_table(star_id,:rrmscdpp10p5)*sqrt(10.5), star_table(star_id,:rrmscdpp12p0)*sqrt(12.), star_table(star_id,:rrmscdpp12p5)*sqrt(12.5), star_table(star_id,:rrmscdpp15p0)*sqrt(15.)]
  contam = 0.0 # rand(LogNormal(1.0e-3,1.0))      # TODO SCI: Come up with better description of Kepler targets, maybe draw from real contaminations
  data_span = star_table(star_id, :dataspan)
  duty_cycle = star_table(star_id, :dutycycle)
  wf_id = 0 # TODO: Set to index containing window function information from stellar table
  # ch = rand(DiscreteUniform(1,84))
  ps = generate_planetary_system(star, sim_param)  
    #return KeplerTarget([ps],fill(cdpp,num_cdpp_timescales,num_quarters),contam,data_span,duty_cycle) #,ch )
  return KeplerTarget([ps],repeat(cdpp_arr, outer=[1,1]),contam,data_span,duty_cycle,wf_id)
end

function generate_kepler_target_simple(sim_param::SimParam)   
  const  generate_star = get_function(sim_param,"generate_star")
  const  generate_planetary_system = get_function(sim_param,"generate_planetary_system")
  const star::StarAbstract = generate_star(sim_param)

  const mean_log_cdpp = 4.9759601617565465   # mean frmo star table
  const stddev_log_cdpp = 0.6704860437536709 # std dev from star table
  rrmscdpp_5hr = exp(mean_log_cdpp+stddev_log_cdpp*randn())
  cdpp = 1.0e-6 * rrmscdpp_5hr * sqrt(5.0/24.0 / LC_duration )
  contam = 0.0 # rand(LogNormal(1.0e-3,1.0))   # TODO SCI: Come up with better description of Kepler targets, maybe draw from real contaminations
  #ch = rand(DiscreteUniform(1,84))
  ps = generate_planetary_system(star, sim_param)  
  return KeplerTarget([ps],fill(cdpp,num_cdpp_timescales,num_quarters),contam,mission_data_span,mission_duty_cycle,0) #,ch)
end

function test_target(sim_param::SimParam)
  generate_kepler_target_simple(sim_param)
  StellarTable.setup_star_table(sim_param)
  generate_kepler_target_from_table(sim_param)
end

