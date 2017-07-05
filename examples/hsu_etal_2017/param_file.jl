using ExoplanetsSysSim

#= 
    p_lim_arr_num = [0.5, 1.25, 2.5, 5., 10., 20., 40., 80., 160., 320.]
    r_lim_arr_num = [1., 1.25, 1.5, 1.75, 2.]
    rate_tab_init = [0.027  0.149 0.354 1.37  2.56  1.96 2.85  4.77 0.0001 ;
                     0.02   0.063 0.432 0.9   1.72  1.46 2.4   5.42 3.63 ;
                     0.0074 0.032 0.389 0.545 1.32  1.65 1.96  1.69 1.17 ;
                     0.0001 0.017 0.156 0.542 0.777 1.35 0.965 1.22 2.49]*0.01
=#

function set_test_param(sim_param_closure::SimParam)
    add_param_fixed(sim_param_closure,"num_targets_sim_pass_one",50000)  # For faster simulated catalogs

    p_lim_arr_num = [10., 20.]
    r_lim_arr_num = [1.0, 1.25]
    rate_tab_init = reshape([2.5]*0.01, (1,1))
    add_param_fixed(sim_param_closure, "p_lim_arr", p_lim_arr_num)
    add_param_fixed(sim_param_closure, "r_lim_arr", r_lim_arr_num*ExoplanetsSysSim.earth_radius)
    add_param_active(sim_param_closure, "obs_par", rate_tab_init)
    return sim_param_closure
end