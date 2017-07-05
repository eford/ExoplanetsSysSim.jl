## ExoplanetsSysSim/src/transit_observations_unused.jl
## (c) 2015 Eric B. Ford

abstract EphemerisAbstract                     # Abstract, so we can generalize to allow for TTVs.  # TODO OPT: Make sure not a big speed penalty
immutable EphemerisLinear <: EphemerisAbstract  # Currently unused
  P::Float64
  t0::Float64
end

immutable TransitShape                       # Currently Unused
  depth::Float64        # Fractional
  duration::Float64     # days               #  QUERY: What definition of duration here?  Since we're going to use expressions from Price & Rogers for measurement uncertainties, let's treat this as the Full-width, half-max-duration  (until further notice)
  b::Float64            # in stellar radii
end

immutable TransitParameter{EphemT}           # Currently Unusued
  ephem::EphemT
  shape::TransitShape
end



if false  # Warning only works if TransitPlanetObs is mutable.  For testing purposes only
function set!(obs::TransitPlanetObs, P::Float64, t0::Float64, d::Float64, D::Float64)
  obs.period = P
  obs.t0 = t0
  obs.depth = d
  obs.duration = D
end
end

