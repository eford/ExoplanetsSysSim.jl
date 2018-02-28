## ExoplanetsSysSim/examples/hsu_etal_2018/invdet_calc.jl
## (c) 2018 Danley C. Hsu

using ExoplanetsSysSim
include(joinpath(Pkg.dir(),"ExoplanetsSysSim","examples","hsu_etal_2018", "christiansen_func.jl"))

global sim_param_closure = setup_sim_param_christiansen()
sim_param_closure = set_test_param(sim_param_closure)

@time inv_det_prob(sim_param_closure)
