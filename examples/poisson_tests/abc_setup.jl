## ExoplanetsSysSim/examples/hsu_etal_2018/abc_setup.jl
## (c) 2018 Eric B. Ford & Danley C. Hsu
# Collection of functions which specific ABC simulation parameters

module EvalSysSimModel
export setup, get_param_vector, get_ss_obs
export gen_data, calc_summary_stats, calc_distance, is_valid, is_valid_single_rp, normalize_dirch
using ExoplanetsSysSim
using ABC
include(joinpath(Pkg.dir(),"ExoplanetsSysSim","examples","poisson_tests", "christiansen_func.jl"))

sim_param_closure = SimParam()
summary_stat_ref_closure =  CatalogSummaryStatistics()

function is_valid(param_vector::Vector{Float64})
    global sim_param_closure
    update_sim_param_from_vector!(param_vector,sim_param_closure)
    const rate_tab::Array{Float64,2} = get_any(sim_param_closure, "obs_par", Array{Float64,2})
    limitP::Array{Float64,1} = get_any(EvalSysSimModel.sim_param_closure, "p_lim_arr", Array{Float64,1})
    #const lambda = sum_kbn(rate_tab)
    if any(x -> x < 0., rate_tab) || any([floor(3*log(limitP[i+1]/limitP[i])/log(2)) for i in 1:length(limitP)-1] .< rate_tab[1,:])
    #if any(x -> x < 0., rate_tab) || any([floor(3*log(limitP[i+1]/limitP[i])/log(2)) for i in 1:length(limitP)-1] .< sum(rate_tab, 1))
        return false
    end
    return true
end

function is_valid_single_rp(param_vector::Vector{Float64})
    global sim_param_closure
    update_sim_param_from_vector!(param_vector,sim_param_closure)
    const rate_tab::Array{Float64,2} = get_any(sim_param_closure, "obs_par", Array{Float64,2})
    limitP::Array{Float64,1} = get_any(EvalSysSimModel.sim_param_closure, "p_lim_arr", Array{Float64,1})
    #const lambda = sum_kbn(rate_tab)
    if any(x -> x < 0., rate_tab) || any([floor(3*log(limitP[i+1]/limitP[i])/log(2)) for i in 1:length(limitP)-1] .< reshape(rate_tab, length(limitP)-1))
        return false
    end
    return true
end

function normalize_dirch(param_vector::Vector{Float64})
    global sim_param_closure
    const p_dim = length(get_any(sim_param_closure, "p_lim_arr", Array{Float64,1}))-1
    const r_dim = length(get_any(sim_param_closure, "r_lim_arr", Array{Float64,1}))-1

    for i in 1:p_dim
        param_vector[((i-1)*(r_dim+1)+2):((i-1)*(r_dim+1)+(r_dim+1))] ./= sum(param_vector[((i-1)*(r_dim+1)+2):((i-1)*(r_dim+1)+(r_dim+1))])
    end

    update_sim_param_from_vector!(param_vector,sim_param_closure)
    return param_vector
end

function gen_data(param_vector::Vector{Float64})
    global sim_param_closure
    update_sim_param_from_vector!(param_vector,sim_param_closure)
    cat_phys = generate_kepler_physical_catalog(sim_param_closure)
    #cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys, sim_param_closure)
    #cat_obs = ExoplanetsSysSim.observe_kepler_targets_single_obs(cat_phys_cut, sim_param_closure)
    cat_obs = ExoplanetsSysSim.observe_kepler_targets_sky_avg(cat_phys, sim_param_closure)
    return cat_obs
end

# TODO OPT: Eventually, could adapt ABC.jl to use distance from first pass to decide if should compute additional summary statistics
function calc_summary_stats(cat::KeplerObsCatalog)
    global sim_param_closure
    sum_stat = calc_summary_stats_obs_binned_rates(cat, sim_param_closure, obs_skyavg = true)
    return sum_stat
end

function calc_distance(sum_stat_obs::CatalogSummaryStatistics,sum_stat_sim::CatalogSummaryStatistics, n::Integer = 0)
    global sim_param_closure
    dist1 = calc_distance_vector_binned(sum_stat_obs,sum_stat_sim, 1, sim_param_closure)
    num_available = length(dist1)
    num_to_use = n>0 ? min(n,num_available) : num_available
    return calc_scalar_distance(dist1[1:num_to_use])
end

function setup()
    global sim_param_closure = setup_sim_param_christiansen()
    sim_param_closure = set_test_param(sim_param_closure)
    
    ### Use simulated planet candidate catalog data
    # df_star = setup_star_table_christiansen(sim_param_closure)
    # println("# Finished reading in stellar data")
    # add_param_fixed(sim_param_closure,"num_kepler_targets",1000000)  # For "observed" catalog
    # cat_obs = simulated_read_kepler_observations(sim_param_closure)
    # println("# Finished setting up simulated true catalog")
    ###
    
    ### Use real planet candidate catalog data
    df_koi,usable_koi = read_koi_catalog(sim_param_closure)
    println("# Finished reading in KOI data")  
    df_star = setup_star_table_christiansen(sim_param_closure)
    println("# Finished reading in stellar data")
    cat_obs = setup_actual_planet_candidate_catalog(df_star, df_koi, usable_koi, sim_param_closure)
    println("# Finished setting up true catalog")
    ###
    
    global summary_stat_ref_closure = calc_summary_stats_obs_binned_rates(cat_obs,sim_param_closure, trueobs_cat = true)
end

get_param_vector() = make_vector_of_sim_param(sim_param_closure)
get_ss_obs() = summary_stat_ref_closure

function set_simparam_ss(sim_param::ExoplanetsSysSim.SimParam, ss_true::ExoplanetsSysSim.CatalogSummaryStatistics)
    global sim_param_closure = sim_param
    global summary_stat_ref_closure = ss_true
end

end  # module EvalSysSimModel

include(joinpath(Pkg.dir("ABC"),"src/composite.jl"))

module SysSimABC
export setup_abc, run_abc, run_abc_largegen, setup_abc_p2
import ABC
import Distributions
using CompositeDistributions
using Compat
import ExoplanetsSysSim
import EvalSysSimModel
include(joinpath(Pkg.dir(),"ExoplanetsSysSim","examples","poisson_tests", "christiansen_func.jl"))
include(joinpath(Pkg.dir(),"ExoplanetsSysSim","examples","poisson_tests", "beta_proposal.jl"))

function setup_abc(num_dist::Integer = 0)
    EvalSysSimModel.setup()
    limitP::Array{Float64,1} = get_any(EvalSysSimModel.sim_param_closure, "p_lim_arr", Array{Float64,1})
    const r_dim = length(get_any(EvalSysSimModel.sim_param_closure, "r_lim_arr", Array{Float64,1}))-1
    prior_arr = ContinuousDistribution[]
    ss_obs_table = EvalSysSimModel.get_ss_obs().stat["planets table"]

    for i in 1:(length(limitP)-1)
        max_in_col = floor(3*log(limitP[i+1]/limitP[i])/log(2))
        lambda_col = Uniform(0.0, max_in_col)
        if r_dim > 1
            dirch = Dirichlet(ones(r_dim))
            prior_arr = vcat(prior_arr, [lambda_col, dirch])
        else
            prior_arr = vcat(prior_arr, lambda_col)
        end
        # for j in 1:r_dim
        #     if j < 4 && limitP[i] > 20.0
        #         push!(prior_arr, Distributions.Uniform(0., max_in_col))
        #     else
        #         push!(prior_arr, Distributions.Uniform(0., 0.2))
        #     end
        # end
    end
    param_prior = CompositeDist(prior_arr)
    in_parallel = nworkers() > 1 ? true : false

    calc_distance_ltd(sum_stat_obs::ExoplanetsSysSim.CatalogSummaryStatistics,sum_stat_sim::ExoplanetsSysSim.CatalogSummaryStatistics) = EvalSysSimModel.calc_distance(sum_stat_obs,sum_stat_sim,num_dist)

    if r_dim > 1
        global abc_plan = ABC.abc_pmc_plan_type(EvalSysSimModel.gen_data,EvalSysSimModel.calc_summary_stats, calc_distance_ltd, param_prior, make_proposal_dist=make_proposal_dist_multidim_beta, is_valid=EvalSysSimModel.is_valid, normalize=EvalSysSimModel.normalize_dirch, num_part=100, num_max_attempt=200, num_max_times=15, epsilon_init=9.9e99, target_epsilon=1.0e-100, in_parallel=in_parallel, adaptive_quantiles = false, epsilon_reduction_factor=0.9, tau_factor=1.5);
    else
        global abc_plan = ABC.abc_pmc_plan_type(EvalSysSimModel.gen_data,EvalSysSimModel.calc_summary_stats, calc_distance_ltd, param_prior, make_proposal_dist=make_proposal_dist_multidim_beta_single_rp, is_valid=EvalSysSimModel.is_valid_single_rp, num_part=100, num_max_attempt=200, num_max_times=15, epsilon_init=9.9e99, target_epsilon=1.0e-100, in_parallel=in_parallel, adaptive_quantiles = false, epsilon_reduction_factor=0.9, tau_factor=1.5);
    end
end

function setup_abc_p2(abc_plan::ABC.abc_pmc_plan_type)
    abc_plan.tau_factor = 2.0
    abc_plan.num_max_times = 200
    ExoplanetsSysSim.add_param_fixed(EvalSysSimModel.sim_param_closure,"num_targets_sim_pass_one",10000)
    return abc_plan
end
    

function run_abc_largegen(pop::ABC.abc_population_type, ss_true::ExoplanetsSysSim.CatalogSummaryStatistics, epshist_targ::Float64; npart::Integer = 1000, num_dist::Integer = 0)
    # sim_param_closure = setup_sim_param_christiansen()
    # sim_param_closure = set_test_param(sim_param_closure)
    # setup_star_table_christiansen(sim_param_closure)
    # EvalSysSimModel.set_simparam_ss(sim_param_closure, ss_true)	
    
    limitP::Array{Float64,1} = get_any(EvalSysSimModel.sim_param_closure, "p_lim_arr", Array{Float64,1})
    const r_dim = length(get_any(EvalSysSimModel.sim_param_closure, "r_lim_arr", Array{Float64,1}))-1
    prior_arr = ContinuousDistribution[]
    for i in 1:(length(limitP)-1)
        max_in_col = floor(3*log(limitP[i+1]/limitP[i])/log(2))
        lambda_col = Uniform(0.0, max_in_col)
        if r_dim > 1
            dirch = Dirichlet(ones(r_dim))
            prior_arr = vcat(prior_arr, [lambda_col, dirch])
        else
            prior_arr = vcat(prior_arr, lambda_col)
        end
        # for j in 1:r_dim
        #     if j < 4 && limitP[i] > 20.0
        #         push!(prior_arr, Distributions.Uniform(0., max_in_col))
        #     else
        #         push!(prior_arr, Distributions.Uniform(0., 0.2))
        #     end
        # end
    end
    param_prior = CompositeDist(prior_arr)
    in_parallel = nworkers() > 1 ? true : false
    
    calc_distance_ltd(sum_stat_obs::ExoplanetsSysSim.CatalogSummaryStatistics,sum_stat_sim::ExoplanetsSysSim.CatalogSummaryStatistics) = EvalSysSimModel.calc_distance(sum_stat_obs,sum_stat_sim,num_dist)

    if r_dim > 1
        global abc_plan = ABC.abc_pmc_plan_type(EvalSysSimModel.gen_data,EvalSysSimModel.calc_summary_stats, calc_distance_ltd, param_prior, make_proposal_dist=make_proposal_dist_multidim_beta, is_valid=EvalSysSimModel.is_valid, normalize=EvalSysSimModel.normalize_dirch, num_part=npart, num_max_attempt=50, num_max_times=1, epsilon_init=9.9e99, target_epsilon=1.0e-100, in_parallel=in_parallel, tau_factor=2.0);
    else
        global abc_plan = ABC.abc_pmc_plan_type(EvalSysSimModel.gen_data,EvalSysSimModel.calc_summary_stats, calc_distance_ltd, param_prior, make_proposal_dist=make_proposal_dist_multidim_beta_single_rp, is_valid=EvalSysSimModel.is_valid_single_rp, num_part=npart, num_max_attempt=50, num_max_times=1, epsilon_init=9.9e99, target_epsilon=1.0e-100, in_parallel=in_parallel, tau_factor=2.0);
    end

    println("# run_abc_largegen: ",EvalSysSimModel.sim_param_closure)
    sampler_largegen = abc_plan.make_proposal_dist(pop, abc_plan.tau_factor)
    theta_largegen = Array{Float64}(size(pop.theta, 1), npart)
    weight_largegen = Array{Float64}(npart)
    for i in 1:npart
        theta_val, dist_largegen, attempts_largegen = ABC.generate_theta(abc_plan, sampler_largegen, ss_true, epshist_targ)
        theta_largegen[:,i] = theta_val  
        prior_logpdf = Distributions.logpdf(abc_plan.prior,theta_val)
        sampler_logpdf = Distributions.logpdf(sampler_largegen, theta_val)
        weight_largegen[i] = exp(prior_logpdf-sampler_logpdf)
    end
    return theta_largegen, weight_largegen
end

function run_abc(abc_plan::ABC.abc_pmc_plan_type)
    #global sim_param_closure
    println("# run_abc: ",EvalSysSimModel.sim_param_closure)
    ss_true = EvalSysSimModel.get_ss_obs()
    #println("True catalog SS: ", ss_true)
    pop_out = ABC.run_abc(abc_plan,ss_true;verbose=true)
end

function run_abc(abc_plan::ABC.abc_pmc_plan_type, pop::ABC.abc_population_type)
    #global sim_param_closure
    dist_threshold = maximum(pop.dist)
    EvalSysSimModel.add_param_fixed(EvalSysSimModel.sim_param_closure,"minimum ABC dist skip pass 2",dist_threshold)
    println("# run_abc: ",EvalSysSimModel.sim_param_closure)
    ss_true = EvalSysSimModel.get_ss_obs()
    pop_out = ABC.run_abc(abc_plan,ss_true,pop;verbose=true)
end

end # module SysSimABC
