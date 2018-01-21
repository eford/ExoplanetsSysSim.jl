using ExoplanetsSysSim

function set_test_param(sim_param_closure::SimParam)
    add_param_fixed(sim_param_closure,"num_targets_sim_pass_one",150969)  # For faster simulated catalogs

    p_lim_arr_num = [237., 320.]
    r_lim_arr_num = [1.0, 1.5]
    rate_tab_init = reshape([2.5]*0.01, (1,1))
    add_param_fixed(sim_param_closure, "p_lim_arr", p_lim_arr_num)
    add_param_fixed(sim_param_closure, "r_lim_arr", r_lim_arr_num*ExoplanetsSysSim.earth_radius)
    add_param_active(sim_param_closure, "obs_par", rate_tab_init)
    return sim_param_closure
end
