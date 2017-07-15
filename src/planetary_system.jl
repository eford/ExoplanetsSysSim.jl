#using Distributions
#include("constants.jl")
#include("orbit.jl")
#include("planet.jl")

if !isdefined(:PlanetarySystemAbstract)
  @compat abstract type PlanetarySystemAbstract end

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
PlanetarySystem{StarT<:StarAbstract}(s::StarT) = PlanetarySystem(s,Array{Planet}(0),Array{Orbit}(0))  # Constructor for a Planetary System with no planets

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
<<<<<<< HEAD
  const mr_max_mass::Float64 = get_real(sim_param,"mr_max_mass")
  m = mr_const*earth_mass*(r/earth_radius)^mr_power_index
  if m > mr_max_mass
     m = mr_max_mass
  end
=======
  m = mr_const*earth_mass*(r/earth_radius)^mr_power_index
>>>>>>> origin/master
  return m
end

function generate_num_planets_poisson(lambda::Real, max_planets_in_sys::Integer; min_planets_in_sys::Integer = 0)
  if lambda < min_planets_in_sys*1e-3
      return min_planets_in_sys
  end
  bug_fixed = false # TODO: true case should work, but Danley found bug in Distributions package.  Revert once fixed for speed.
  local n
  if bug_fixed
     d = Distributions.Truncated(Distributions.Poisson(lambda),min_planets_in_sys,max_planets_in_sys)
     n = rand(d)
  else
     if min_planets_in_sys == 0 
        min_planets_in_sys = -1
     end 
     d = Distributions.Truncated(Distributions.Poisson(lambda),min_planets_in_sys,max_planets_in_sys)
     n = rand(d)
     #=
     n = -1
     while !(min_planets_in_sys<=n<=max_planets_in_sys)
        n = rand(Distributions.Poisson(lambda))
     end 
     =#
  end
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

<<<<<<< HEAD
function generate_periods_power_law(s::Star, sim_param::SimParam; num_pl::Integer = 1) 
    const power_law_P::Float64 = get_real(sim_param,"power_law_P")
    const min_period::Float64 = get_real(sim_param,"min_period")
    const max_period::Float64 = get_real(sim_param,"max_period")
    Plist = power_law_P!=-1.0 ? draw_power_law(power_law_P,min_period,max_period, num_pl) : exp(log(min_period).+rand()*log(max_period/min_period))
    return Plist
end

function generate_sizes_power_law(s::Star, sim_param::SimParam; num_pl::Integer = 1) 
    const power_law_r::Float64 = get_real(sim_param,"power_law_r")
    const min_radius::Float64 = get_real(sim_param,"min_radius")
    const max_radius::Float64 = get_real(sim_param,"max_radius")
    Rlist = power_law_r!=-1.0 ? draw_power_law(power_law_r,min_radius,max_radius, num_pl) : exp(log(min_radius).+rand()*log(max_radius/min_radius))
    return Rlist
end

function generate_period_and_sizes_power_law(s::Star, sim_param::SimParam; num_pl::Integer = 1) 
    return generate_periods_power_law(s, sim_param, num_pl=num_pl), generate_sizes_power_law(s, sim_param, num_pl=num_pl)
end


=======
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

>>>>>>> origin/master
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
     pl = Array{Planet}(num_pl)
     orbit = Array{Orbit}(num_pl)
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
     pl = Array{Planet}(num_pl)
     orbit = Array{Orbit}(num_pl)
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
     pl = Array{Planet}(num_pl)
     orbit = Array{Orbit}(num_pl)
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

<<<<<<< HEAD
# NEW STUFF FOR Matthias BELOW

function calc_hill_sphere(a::Float64, mu::Float64)
  a*(mu/3)^(1//3)
end

function calc_mutual_hill_radii{StarT<:StarAbstract}(ps::PlanetarySystem{StarT},pl1::Int64, pl2::Int64)
  mu = (ps.planet[pl1].mass+ps.planet[pl2].mass)/ps.star.mass
  a = 0.5*(ps.orbit[pl1].a+ps.orbit[pl2].a)
  calc_hill_sphere(a,mu)
end

function test_stability_circular(P::Vector{Float64},mass::Vector{Float64},star_mass::Float64)
   @assert length(P) == length(mass)
   const min_num_mutual_hill_radii = get_real(sim_param,"num_mutual_hill_radii")
   found_instability = false
   order = sortperm(P)
   a2 = semimajor_axis(P[order[1]],star_mass)
   for pl in 1:(length(P)-1)
       a1 = a2   # semimajor_axis(P[order[pl]],star_mass)
       a2 = semimajor_axis(P[order[pl+1]],star_mass)
       a = 0.5*(a1+a2)
       mu = (mass[order[pl]]+mass[order[pl+1]])/star_mass
       mutual_hill_radius = calc_hill_sphere(a,mu)
       if ! ( a2-a1  >= min_num_mutual_hill_radii*mutual_hill_radius )
          found_instability = true
          break
       end
   end # loop over neighboring planet pairs within cluster
   return !found_instability
end

function test_stability(P::Vector{Float64},mass::Vector{Float64},star_mass::Float64; ecc::Vector{Float64} = zeros(length(P)) )
   @assert length(P) == length(mass) == length(ecc)
   const min_num_mutual_hill_radii = get_real(sim_param,"num_mutual_hill_radii")
   found_instability = false
   order = sortperm(P)
   a2 = semimajor_axis(P[order[1]],star_mass)
   for pl in 1:(length(P)-1)
       a1 = a2   # semimajor_axis(P[order[pl]],star_mass)
       a2 = semimajor_axis(P[order[pl+1]],star_mass)
       a = 0.5*(a1+a2)
       mu = (mass[order[pl]]+mass[order[pl+1]])/star_mass
       mutual_hill_radius = calc_hill_sphere(a,mu)
       e1 = ecc[order[pl]]
       e2 = ecc[order[pl+1]]
       if ! ( a2*(1-e2)-a1*(1+e1)  >= min_num_mutual_hill_radii*mutual_hill_radius )
          found_instability = true
          break
       end
   end # loop over neighboring planet pairs within cluster
   return !found_instability
end

function calc_if_near_resonance(P::Vector{Float64})
   @assert issorted(P)   # TODO: OPT: Could remove once know it's safe
   const resonance_width = 0.05   # TODO: FEATURE Make a model parameter?
   const resonance_width_factor = 1+resonance_width
   const period_ratios_to_check = [ 2.0, 1.5, 4/3, 5/4 ]
   result = false(length(P))
   if length(P) >= 2 
      for i in 1:(length(P)-1)
          for period_ratio in period_ratios_to_check
              if P[i]*period_ratio <= P[i+1] <= P[i]*period_ratio*resonance_width_factor
                 result[i] = true
                 result[i+1] = true
                 break 
              end # if near resonance
          end # period_ratio
      end # planets
   end # at least two planets
   return result
end

function generate_planet_periods_sizes_masses_in_cluster( star::StarAbstract, sim_param::SimParam; n::Int64 = 1 )  # TODO: IMPORTANT: Make this function work and test before using for science
   @assert n>=1
   const generate_planet_mass_from_radius = get_function(sim_param,"generate_planet_mass_from_radius")

   if n==1
      R = generate_sizes_power_law(star,sim_param)
      mass = map(generate_planet_mass_from_radius,R)
      P = [1.0] # generate_periods_power_law(star,sim_param)
      return P, R, mass
   end

   # If reach here, then at least 2 planets in cluster

   mean_R = generate_sizes_power_law(star,sim_param)
   const sigma_log_radius_in_cluster = get_real(sim_param,"sigma_log_radius_in_cluster")
   Rdist = LogNormal(log(mean_R),sigma_log_radius_in_cluster)
   R = rand(Rdist,n)

   mass = map(generate_planet_mass_from_radius,R)

   const sigma_logperiod_per_pl_in_cluster = get_real(sim_param,"sigma_logperiod_per_pl_in_cluster")
   log_mean_P = 0.0 # log(generate_periods_power_law(star,sim_param))
   # Note: Currently, drawing all periods within a cluster at once and either keeping or rejecting the whole cluster
   #       Should we instead draw periods one at a time?
   Pdist = LogNormal(log_mean_P,sigma_logperiod_per_pl_in_cluster*n)
   local P
   found_good_periods = false
   attempts = 0
   max_attempts = 100  # Note: Should this be a parameter?
   while !found_good_periods && attempts<max_attempts
      attempts += 1
      P = rand(Pdist,n)
      if test_stability_circular(P,mass,star.mass)
        found_good_periods = true
     end
   end # while trying to draw periods 
   if !found_good_periods
      println("# Warning: Didn't find a good set of periods, sizes and masses for one cluster.")
      return fill(NaN,n), R, mass  # Return NaN's for periods to indicate failed
   end
   return P, R, mass    # Note can also return earlier if only one planet in cluster or if fail to generate a good set of values
end



# This version generates clustered planetary systems, adaptation of Matthias He's python code
function generate_planetary_system_clustered(star::StarAbstract, sim_param::SimParam; verbose::Bool = false)  # TODO: Make this function work and test before using for science
  # load functions to use for drawing parameters
  const generate_num_clusters = get_function(sim_param,"generate_num_clusters")
  const generate_num_planets_in_cluster = get_function(sim_param,"generate_num_planets_in_cluster")
  
  # Generate a set of periods, planet radii, and planet masses.
  attempt_system = 0
  max_attempts_system = 100
  valid_system = false
  while !valid_system && attempt_system <= max_attempts_system

     # First, generate number of clusters (to attempt) and planets (to attempt) in each cluster
     num_clusters = generate_num_clusters(star,sim_param)::Int64
     num_pl_in_cluster = map(x->generate_num_planets_in_cluster(star, sim_param)::Int64, 1:num_clusters)
     num_pl = sum(num_pl_in_cluster)
 
     if( num_pl==0 )
       return PlanetarySystem(star)
     end
     Plist = Array{Float64}(num_pl)
     Rlist = Array{Float64}(num_pl)
     masslist = Array{Float64}(num_pl)
     @assert num_pl_in_cluster[1] >= 1
     pl_start = 1
     pl_stop = 0
     for c in 1:num_clusters
         pl_stop  += num_pl_in_cluster[c]
         Plist_tmp, Rlist[pl_start:pl_stop], masslist[pl_start:pl_stop] = generate_planet_periods_sizes_masses_in_cluster(star, sim_param, n = num_pl_in_cluster[c])
         valid_cluster = !any(isnan.(Plist))
         valid_period_scale = false
         attempt_period_scale = 0
         max_attempts_period_scale = 100
         while !valid_period_scale && attempt_period_scale<=max_attempts_period_scale && valid_cluster
             attempt_period_scale += 1
             period_scale = generate_periods_power_law(star,sim_param)
             Plist[pl_start:pl_stop] = Plist_tmp * period_scale
             if test_stability_circular(Plist[1:pl_stop],masslist[1:pl_stop],star.mass)  # Note: Should we include eccentricities in this test?
                valid_period_scale = true
             end
         end  # while !valid_period_scale...
         if !valid_period_scale
            Plist[pl_star:pl_stop] = NaN
         end
         pl_start += num_pl_in_cluster[c]
     end # for c in 1:num_clusters
     if any(isnan.(Plist))  # If any loop failed to generate valid planets, it should set a NaN in the period list
        keep = !(isnan.(Plist))  # Currently, keeping clusters that could be fit, rather than throwing out entirely and starting from scratch.  Is this a good idea?  Matthias tried the other approach in his python code.  
        num_pl = sum(keep)
        Plist = Plist[keep]
        Rlist = Rlist[keep]
        masslist = masslist[keep]
     end
     valid_system = true
     #= Note: This would be for drawing each cluster separately and then accepting or rejecting the whole lot.
             By testing for stability before adding each cluster, this last test should be unnecessary.
     if test_stability(Plist,masslist,star.mass)
        valid_system = true
     end
     =#
     attempt_system += 1
  end # while !valid_system...

  # Now assign orbits with given periods, sizes, and masses.
  const generate_e_omega = get_function(sim_param,"generate_e_omega")
  const sigma_ecc::Float64 = haskey(sim_param,"sigma_hk") ? get_real(sim_param,"sigma_hk") : 0.0
  const sigma_incl = deg2rad(get_real(sim_param,"sigma_incl"))
  const sigma_incl_near_mmr = deg2rad(get_real(sim_param,"sigma_incl_near_mmr"))
    pl = Array{Planet}(num_pl)
    orbit = Array{Orbit}(num_pl)
    incl_sys = acos(rand())
    idx = sortperm(Plist)       # TODO OPT: Check to see if sorting is significant time sink.  If so, could reduce redundant sortperm
    is_near_resonance = calc_if_near_resonance(Plist[idx])
    for i in 1:num_pl
      # if verbose   println("i=",i," idx=",idx," Plist=",Plist[idx] );     end
      if haskey(sim_param,"sigma_hk_one") && haskey(sim_param,"sigma_hk_multi")
         sigma_ecc = num_pl == 1 ? get_real(sim_param,"sigma_hk_one") : get_real(sim_param,"sigma_hk_multi")
      end
      (ecc,  omega) = generate_e_omega(sim_param)::Tuple{Float64,Float64}  # WARNING: Not testing for stability after eccentricites drawn.
      sigma_incl_use = is_near_resonance[i] ? sigma_incl_near_mmr : sigma_incl
      incl_mut = sigma_incl_use*sqrt(randn()^2+randn()^2) # rand(Distributions.Rayleigh(sigma_incl)) 
      asc_node = 2pi*rand()
      mean_anom = 2pi*rand()
      incl =  incl_mut!=zero(incl_mut) ? acos( cos(incl_sys)*cos(incl_mut) + sin(incl_sys)*sin(incl_mut)*cos(asc_node) ) : incl_sys
      orbit[idx[i]] = Orbit(Plist[idx[i]],ecc,incl,omega,asc_node,mean_anom)
      pl[i] = Planet( Rlist[idx[i]],  masslist[idx[i]] )
    end # for i in 1:num_pl

  return PlanetarySystem(star,pl,orbit)
end



=======
>>>>>>> origin/master
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
<<<<<<< HEAD

=======
>>>>>>> origin/master
