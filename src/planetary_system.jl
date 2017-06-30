#using Distributions
#include("constants.jl")
#include("orbit.jl")
#include("planet.jl")

if !isdefined(:PlanetarySystemAbstract)
  abstract PlanetarySystemAbstract

  immutable PlanetarySystem{StarT<:StarAbstract} <: PlanetarySystemAbstract
    star::StarT                    
    planet::Vector{Planet}         # TODO OPT: If want to optimize, consider using  something like ImmmutableArrays?
    orbit::Vector{Orbit}           # TODO OPT: If want to optimize, consider using  something like ImmmutableArrays?
    # TODO DETAIL: Setup inner constructor to enforce equal number of planets & orbits
    #PlanetarySystem{StarT}(s::StarAbstract, p::Vector{Planet}, o::Vector{Orbit})
    #  @assert(length(p)==length(o)) # else error(string("Number of planets must match number of orbits: Np= ",length(p)," No= ",length(o)))
    #  new(s,p,o)
    #end
  end
  #typealias PlanetarySystemSingleStar  PlanetarySystem{SingleStar}
  PlanetarySystemSingleStar = PlanetarySystem{SingleStar}

end
PlanetarySystem{StarT<:StarAbstract}(s::StarT) = PlanetarySystem(s,Array(Planet,0),Array(Orbit,0))  # Constructor for a Planetary System with no planets

PlanetarySystem{StarT<:StarAbstract}(s::StarT, p::Planet, o::Orbit) = PlanetarySystem(s,[p],[o])  # Constructor for a single Planet System

function PlanetarySystem{StarT}(ps::PlanetarySystem{StarT}, keep::Vector{Int64}) # Why doesn't this work?
   PlanetarySystem{StarT}(ps.star,ps.planet[keep],ps.orbit[keep])
end

#function PlanetarySystemSingleStar(ps::PlanetarySystemSingleStar, keep::Vector{Int64})
#function PlanetarySystem{StarT<:StarAbstract}(ps::PlanetarySystem{StarT}, keep::Vector{Int64})
#   PlanetarySystem(ps.star,ps.planet[keep],ps.orbit[keep])
#end

flux(ps::PlanetarySystem{StarAbstract}) = flux(ps.star)
flux(ps::PlanetarySystem{Star}) = flux(ps.star)
flux(ps::PlanetarySystem{BinaryStar}) = flux(ps.star)
flux(ps::PlanetarySystem{MultipleStar}) = flux(ps.star)

function num_planets(s::PlanetarySystem{Star})
  @assert( length(s.planet) == length(s.orbit) )    # TODO OPT: Deactivate inner assert's like this for speed once tested
  return length(s.planet)
end

function generate_planet_mass_from_radius_powerlaw(r::Float64, s::Star, o::Orbit, sim_param::SimParam) # TODO SCI: IMPORTANT once have stability criteria, replace w/ better M-R relationship
  const mr_power_index::Float64 = get_real(sim_param,"mr_power_index")
  const mr_const::Float64 = get_real(sim_param,"mr_const")
  m = mr_const*earth_mass*(r/earth_radius)^mr_power_index
  return m
end

function generate_num_planets_poisson(lambda::Real, max_planets_in_sys::Integer; min_planets_in_sys::Integer = 0)
  if lambda < min_planets_in_sys*1e-3
      return min_planets_in_sys
  end
  d = Distributions.Truncated(Distributions.Poisson(lambda),min_planets_in_sys,max_planets_in_sys)
  n = rand(d)
  #=
  n = -1
  while !(min_planets_in_sys<=n<=max_tranets_in_sys)
     n = rand(Distributions.Poisson(lambda))
  end 
  =#
 return n
end

function generate_num_planets_poisson(s::Star, sim_param::SimParam)
  const lambda::Float64 = exp(get_real(sim_param,"log_eta_pl"))
  const max_tranets_in_sys::Int64 = get_int(sim_param,"max_tranets_in_sys")
  generate_num_planets_poisson(lambda,max_tranets_in_sys)
end

function generate_period_and_sizes_log_normal(s::Star, sim_param::SimParam; num_pl::Integer = 1)  # TODO SCI:  IMPORTANT: Improve how periods and sizes are drawn
    const mu_log_r::Float64 = get_real(sim_param,"mean_log_planet_radius")
    const sigma_log_r::Float64 = get_real(sim_param,"sigma_log_planet_radius")
    const mu_log_P::Float64 = get_real(sim_param,"mean_log_planet_period")
    const sigma_log_P::Float64 = get_real(sim_param,"sigma_log_planet_period")
    const min_period::Float64 = get_real(sim_param,"min_period")
    const max_period::Float64 = get_real(sim_param,"max_period")
    const min_radius::Float64 = get_real(sim_param,"min_radius")
    const max_radius::Float64 = get_real(sim_param,"max_radius")
    const max_draws::Int64 = 100
 
    if   sigma_log_r <= 0. || sigma_log_P<=0.
     println("# mu_log_r= ", mu_log_r, " sigma_log_r= ", sigma_log_r, " mu_log_P= ", mu_log_P, " sigma_log_P= ", sigma_log_P)
    end
    const rdist = LogNormal(mu_log_r,sigma_log_r)
    const Pdist = LogNormal(mu_log_P,sigma_log_P)
    #Rlist = rand(rdist,num_pl)
    #Plist = rand(Pdist,num_pl)
    #idx_keep = find(i->(min_radius<=Rlist[i]<=max_radius) && (min_period<=Plist[i]<=max_period), 1:num_pl )
    #return Plist[idx_keep], Rlist[idx_keep]  # replaced because want to return exactly num_pl
    Rlist = zeros(num_pl)
    Plist = zeros(num_pl)
    for i in 1:num_pl
      j = 0
      while ! (min_radius<Rlist[i]<=max_radius) && j<max_draws
            Rlist[i] = rand(rdist)
            j+=1
      end
      if j>=max_draws
         println("# Struggled to draw size for: ",mu_log_r, " ", sigma_log_r)
      end
      j = 0
      while ! (min_period<Plist[i]<=max_period) && j<max_draws
            Plist[i] = rand(Pdist)
            j+=1
      end
      if j>=max_draws
         println("# Struggled to draw period for: ",mu_log_P, " ", sigma_log_P)
      end
    end
    return Plist, Rlist
end

function draw_power_law(n::Real, x0::Real, x1::Real, num_pl::Integer)
     ((x1^(n+1) - x0^(n+1))*rand(num_pl) + x0^(n+1)).^(1/(n+1))
end

function generate_period_and_sizes_power_law(s::Star, sim_param::SimParam; num_pl::Integer = 1) 
    const power_law_r::Float64 = get_real(sim_param,"power_law_r")
    const power_law_P::Float64 = get_real(sim_param,"power_law_P")
    const min_period::Float64 = get_real(sim_param,"min_period")
    const max_period::Float64 = get_real(sim_param,"max_period")
    const min_radius::Float64 = get_real(sim_param,"min_radius")
    const max_radius::Float64 = get_real(sim_param,"max_radius")
    Rlist = power_law_r!=-1.0 ? draw_power_law(power_law_r,min_radius,max_radius, num_pl) : exp(log(min_radius).+rand()*log(max_radius/min_radius))
    Plist = power_law_P!=-1.0 ? draw_power_law(power_law_P,min_period,max_period, num_pl) : exp(log(min_period).+rand()*log(max_period/min_period))
    return Plist, Rlist
end

function generate_e_omega_rayleigh(sigma_hk::Float64)
  h = k = 1.0
  while h*h+k*k >= 1.0
    h = sigma_hk*randn()
    k = sigma_hk*randn()
  end
  ecc = sqrt(h*h+k*k)
  w = atan2(k,h)
  return ecc, w
end

function generate_e_omega_rayleigh(sim_param::SimParam)
  sigma_hk::Float64 = get_real(sim_param,"sigma_hk")
  generate_e_omega_rayleigh(sigma_hk)
end

function generate_planetary_system_hardcoded_example(star::StarAbstract, sim_param::SimParam; verbose::Bool = false)
  # in this version we specify fixed functions that are known at compile time, allowing for additional optimizations (~0.6 second faster per Kepler catalog out of ~3.6 sec on my laptop w/ 1 core)
  const  generate_planet_mass_from_radius = generate_planet_mass_from_radius_powerlaw
  const  generate_num_planets = generate_num_planets_poisson
  const  generate_period_and_sizes = generate_period_and_sizes_log_normal

  const  generate_e_omega =  generate_e_omega_rayleigh

  # const  generate_star = get_function(sim_param,"generate_star")
  # const star::StarAbstract = generate_star(sim_param)
  const num_pl::Int64 = generate_num_planets(star, sim_param)

  if( num_pl==0 )
    return PlanetarySystem(star)
  else
     pl = Array(Planet,num_pl)
     orbit = Array(Orbit,num_pl)
     (Plist::Vector{Float64}, Rlist::Vector{Float64}) = generate_period_and_sizes(star, sim_param, num_pl=num_pl)
     idx = sortperm(Plist)                   # TODO OPT: Check to see if sorting is significant time sink.  If so, it might could be deferred

    for i in 1:num_pl
      # if verbose   println("i=",i," idx=",idx," Plist=",Plist[idx] );     end
      P = Plist[idx[i]]
      (ecc::Float64,  omega::Float64) = generate_e_omega(sim_param)
      incl::Float64 = acos(rand())
      orbit[idx[i]] = Orbit(P,ecc,incl,omega,2pi*rand(),2pi*rand())
      # set!(orbit[idx[i]],P,ecc,incl,omega,2pi*rand(),2pi*rand())
      mass::Float64 = generate_planet_mass_from_radius(Rlist[idx[i]], star, orbit[idx[i]], sim_param)
      pl[i] = Planet( Rlist[idx[i]],  mass )
    end
  return PlanetarySystem(star,pl,orbit)
  end
end

function generate_planetary_system_uncorrelated_incl(star::StarAbstract, sim_param::SimParam; verbose::Bool = false)
  # load functions to use for drawing parameters
  const  generate_planet_mass_from_radius = get_function(sim_param,"generate_planet_mass_from_radius")
  const  generate_num_planets = get_function(sim_param,"generate_num_planets")
  const  generate_period_and_sizes = get_function(sim_param,"generate_period_and_sizes")
  const  generate_e_omega = get_function(sim_param,"generate_e_omega")

  # const  generate_star = get_function(sim_param,"generate_star")
  # const star::StarAbstract = generate_star(sim_param)
  const num_pl::Int64 = generate_num_planets(star, sim_param)::Int64
  const sigma_ecc::Float64 = haskey(sim_param,"sigma_hk") ? get_real(sim_param,"sigma_hk") : 0.0

  if( num_pl==0 )
    return PlanetarySystem(star)::PlanetarySystem
  else
     pl = Array(Planet,num_pl)
     orbit = Array(Orbit,num_pl)
     (Plist::Vector{Float64}, Rlist::Vector{Float64}) = generate_period_and_sizes(star, sim_param, num_pl=num_pl)
     idx = sortperm(Plist)                   # TODO OPT: Check to see if sorting is significant time sink.  If so, it might could be deferred

    for i in 1:num_pl
      # if verbose   println("i=",i," idx=",idx," Plist=",Plist[idx] );     end
      P = Plist[idx[i]]
      if haskey(sim_param,"sigma_hk_one") && haskey(sim_param,"sigma_hk_multi")
         sigma_ecc = num_pl == 1 ? get_real(sim_param,"sigma_hk_one") : get_real(sim_param,"sigma_hk_multi")
      end
      (ecc::Float64,  omega::Float64) = generate_e_omega(sim_param)
      incl::Float64 = acos(rand())
      orbit[idx[i]] = Orbit(P,ecc,incl,omega,2pi*rand(),2pi*rand())
      # set!(orbit[idx[i]],P,ecc,incl,omega,2pi*rand(),2pi*rand())
      mass::Float64 = generate_planet_mass_from_radius(Rlist[idx[i]], star, orbit[idx[i]], sim_param)
      pl[i] = Planet( Rlist[idx[i]],  mass )
    end
  return PlanetarySystem(star,pl,orbit)
  end
end

# This version generates more systems roughly near a common plane, but until incorporate CORBITS data, ABC can't match input param
function generate_planetary_system_simple(star::StarAbstract, sim_param::SimParam; verbose::Bool = false)
  # load functions to use for drawing parameters
  const  generate_planet_mass_from_radius = get_function(sim_param,"generate_planet_mass_from_radius")
  const  generate_num_planets = get_function(sim_param,"generate_num_planets")
  #const  generate_num_planets = generate_num_planets_christiansen
  #const  generate_period_and_sizes = generate_period_and_sizes_christiansen
  const  generate_period_and_sizes = get_function(sim_param,"generate_period_and_sizes")
  const  generate_e_omega = get_function(sim_param,"generate_e_omega")
  const  sigma_incl = deg2rad(get_real(sim_param,"sigma_incl"))

  # const  generate_star = get_function(sim_param,"generate_star")
  # const star::StarAbstract = generate_star(sim_param)
  const num_pl = generate_num_planets(star, sim_param)::Int64
  const sigma_ecc::Float64 = haskey(sim_param,"sigma_hk") ? get_real(sim_param,"sigma_hk") : 0.0

  if( num_pl==0 )
    return PlanetarySystem(star)
  else
     pl = Array(Planet,num_pl)
     orbit = Array(Orbit,num_pl)
     (Plist::Vector{Float64}, Rlist::Vector{Float64}) = generate_period_and_sizes(star, sim_param, num_pl=num_pl)
     idx = sortperm(Plist)                   # TODO OPT: Check to see if sorting is significant time sink.  If so, it might could be deferred
     incl_sys = acos(rand())

    for i in 1:num_pl
      # if verbose   println("i=",i," idx=",idx," Plist=",Plist[idx] );     end
      P = Plist[idx[i]]
      if haskey(sim_param,"sigma_hk_one") && haskey(sim_param,"sigma_hk_multi")
         sigma_ecc = num_pl == 1 ? get_real(sim_param,"sigma_hk_one") : get_real(sim_param,"sigma_hk_multi")
      end
      (ecc,  omega) = generate_e_omega(sim_param)::Tuple{Float64,Float64}
      incl_mut = sigma_incl*sqrt(randn()^2+randn()^2) # rand(Distributions.Rayleigh(sigma_incl)) # sigma_incl*randn()
      asc_node = 2pi*rand()
      mean_anom = 2pi*rand()
      #incl = incl_sys + sigma_incl*randn()
      incl =  incl_mut!=zero(incl_mut) ? acos( cos(incl_sys)*cos(incl_mut) + sin(incl_sys)*sin(incl_mut)*cos(asc_node) ) : incl_sys
      orbit[idx[i]] = Orbit(P,ecc,incl,omega,asc_node,mean_anom)
      # set!(orbit[idx[i]], P,ecc,incl,omega,asc_node,mean_anom) # if Orbit were mutable
      mass = generate_planet_mass_from_radius(Rlist[idx[i]], star, orbit[idx[i]], sim_param)::Float64
      pl[i] = Planet( Rlist[idx[i]],  mass )
    end
  return PlanetarySystem(star,pl,orbit)
  end
end

function test_planetary_system_constructors(sim_param::SimParam)
  generate_star = get_function(sim_param,"generate_star")
  star = generate_star(sim_param)
  empty_sys = PlanetarySystem(star)
  earth = Planet(earth_radius,earth_mass)
  earth_orbit = Orbit(365.2425,0.0167,0.5*pi,0.0,0.0,0.0)
  solar_sys = PlanetarySystem(star, earth,earth_orbit)
  m = generate_planet_mass_from_radius_powerlaw(0.02,star,earth_orbit,sim_param)/earth_mass
  generate_planetary_system_simple(star,sim_param,verbose=true)
end
