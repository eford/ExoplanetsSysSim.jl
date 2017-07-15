using ExoplanetsSysSim

function set_test_param(sim_param_closure::SimParam)
    p_lim_arr_num = [0.5, 1.25, 2.5, 5., 10., 20., 40., 80., 160., 320.]
    r_lim_arr_num = [0.5, 0.75, 1., 1.25, 1.5, 1.75, 2., 2.5, 3., 4., 6., 8., 12., 16.]
    rate_tab_init = [0.004 0.031 0.124 0.075 0.162 0.175 0.328 0.503 0.591 ;
                     0.004 0.012 0.082 0.165 0.234 0.407 0.321 1.223 2.381 ;
                     0.001 0.007 0.057 0.103 0.186 0.234 0.764 0.700 1.915 ;
                     0.009 0.012 0.127 0.229 0.483 0.822 0.910 1.463 1.795 ;
                     0.004 0.029 0.136 0.599 1.409 1.889 2.279 3.366 4.987 ;
                     0.007 0.042 0.218 0.932 2.270 2.832 3.500 3.135 3.166 ;
                     0.010 0.059 0.289 1.148 2.527 3.658 3.550 3.616 6.367 ;
                     0.036 0.057 0.380 0.796 1.552 1.181 1.959 2.420 3.505 ;
                     0.047 0.126 0.738 1.089 1.909 1.924 2.254 2.675 2.199 ;
                     0.053 0.126 0.630 1.608 2.405 2.411 2.272 4.370 7.607 ;
                     0.088 0.180 0.524 1.523 2.801 3.055 4.071 3.216 8.540 ;
                     0.100 0.295 1.022 3.156 3.559 7.929 0.0852 2.722 10.4 ;
                     0.039 0.404 2.615 5.095 5.378 1.969 7.989  0.0   0.0  ] * 0.01
    add_param_fixed(sim_param_closure, "p_lim_arr", p_lim_arr_num)
    add_param_fixed(sim_param_closure, "r_lim_arr", r_lim_arr_num*ExoplanetsSysSim.earth_radius)
    add_param_fixed(sim_param_closure, "obs_par", rate_tab_init)

  add_param_fixed(sim_param_closure,"sigma_hk",0.03)
  add_param_fixed(sim_param_closure,"sigma_incl",1.0)   # degrees; 0 = coplanar w/ generate_kepler_target_simple; ignored by generate_planetary_system_uncorrelated_incl
 
    return sim_param_closure
end

