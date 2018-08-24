## ExoplanetsSysSim/src/limb_darkening.jl
## (c) 2018 Eric B. Ford

abstract type LimbDarkeningParamAbstract end
# Note this file is currently not used by SysSim.

struct LimbDarkeningParamLinear <: LimbDarkeningParamAbstract
  coeff::Tuple{Float64}
  # TODO SCI DETAIL: Repalce with sensible limits on LD params
  LimbDarkeningParamLinear(_u1::Real ) = (!( (-2.0<=_u1<=2.0) )) ? error(string("Invalid quadratic limb darkening parameters: ",_u1," & ", _u2)) : new((_u1,))
end

function depth_at_midpoint(radius_ratio::Float64, ld::LimbDarkeningParamLinear)
    c0 = 1-sum(ld.coeff)  
    omega = c0/4+ld.coeff[1]/6
    ksq = 1-radius_ratio^2
    tmp0 = c0/4*ksq
    tmp2 = ld.coeff[1]/6*ksq^(3//2)
    return 1- (tmp0+tmp2)/omega
end

struct LimbDarkeningParamQuadratic <: LimbDarkeningParamAbstract
  coeff::Tuple{Float64,Float64}
  # TODO SCI DETAIL: Repalce with sensible limits on LD params
  LimbDarkeningParamQuadratic(_u1::Real, _u2::Real ) = (!( (-2.0<=_u1<=2.0) && (-2.0<=_u2<=2.0) )) ? error(string("Invalid quadratic limb darkening parameters: ",_u1," & ", _u2)) : new((_u1,_u2))
end

function depth_at_midpoint(radius_ratio::Float64, ld::LimbDarkeningParamQuadratic)
    c0 = 1-sum(ld.coeff) 
    omega = c0/4+(ld.coeff[1]+2*ld.coeff[2])/6-ld.coeff[2]/8 
    ksq = 1-radius_ratio^2
    tmp0 = c0/4*ksq
    tmp2 = (ld.coeff[1]+2*ld.coeff[2])/6*ksq^(3//2)
    tmp4 = -ld.coeff[2]/8*ksq^2
    return 1- (tmp0+tmp2+tmp4)/omega
end

struct LimbDarkeningParam4thOrder <: LimbDarkeningParamAbstract
  coeff::Tuple{Float64,Float64,Float64,Float64}
  # TODO SCI DETAIL: Repalce with sensible limits on LD params
  LimbDarkeningParam4thOrder(_c1::Real, _c2::Real, _c3::Real, _c4::Real ) = (!( (-2.0<=_c1<=2.0) && (-2.0<=_c2<=2.0)  && (-2.0<=_c3<=2.0) && (-2.0<=_c4<=2.0) )) ? error(string("Invalid limb darkening parameters: ",_c1," , ", _c2, ", ",_c3," , ",_c4)) : new((_c1,_c2,_c3,_c4))
end

function depth_at_midpoint(radius_ratio::Float64, ld::LimbDarkeningParam4thOrder)
    c0 = 1-sum(ld.coeff)  
    omega = c0/4+sum(ld.coeff./(5:8))
    ksq = 1-radius_ratio^2
    tmp0 = c0/4*ksq
    tmp1 = ld.coeff[1]/5*ksq^(5//4)
    tmp2 = ld.coeff[2]/6*ksq^(3//2)
    tmp3 = ld.coeff[3]/7*ksq^(7//4)
    tmp4 = ld.coeff[4]/8*ksq^2
    return 1- (tmp0+tmp1+tmp2+tmp3+tmp4)/omega
end

