using ExoplanetsSysSim
include("christiansen_func.jl")

global sim_param_closure = setup_sim_param_christiansen()

p_lim_arr_num = [237, 320.]
r_lim_arr_num = [1., 1.5]

add_param_fixed(sim_param_closure, "p_lim_arr", p_lim_arr_num)
add_param_fixed(sim_param_closure, "r_lim_arr", r_lim_arr_num*ExoplanetsSysSim.earth_radius)

@time inv_det_prob(sim_param_closure)
