using ExoplanetsSysSim
include("christiansen_func.jl")

global sim_param_closure = setup_sim_param_christiansen()

#p_lim_arr_num = [0.5, 1.25, 2.5, 5., 10., 20., 40., 80., 160., 320.]
#r_lim_arr_num = [1., 1.25, 1.5, 1.75, 2.]

p_lim_arr_num = [10., 20.]
r_lim_arr_num = [1., 1.25]

add_param_fixed(sim_param_closure, "p_lim_arr", p_lim_arr_num)
add_param_fixed(sim_param_closure, "r_lim_arr", r_lim_arr_num*ExoplanetsSysSim.earth_radius)

@time inv_det_prob(sim_param_closure)