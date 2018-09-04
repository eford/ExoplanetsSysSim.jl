## ExoplanetsSysSim/examples/hsu_etal_2018/abc_setup.jl
## (c) 2018 Eric B. Ford & Danley C. Hsu
# Collection of functions which specific ABC simulation parameters

module EvalSysSimModel
export setup, get_param_vector, get_ss_obs #, evaluate_model
export gen_data, calc_summary_stats, calc_distance, is_valid, normalize_dirch
export make_proposal_dist_multidim_beta
using ExoplanetsSysSim
using ABC
using SpecialFunctions
using Compat.Statistics
import ABC.CompositeDistributions.CompositeDist
import ABC.TransformedBetaDistributions.LinearTransformedBeta
include(joinpath(Pkg.dir(),"ExoplanetsSysSim","examples","poisson_tests", "christiansen_func.jl"))

sim_param_closure = SimParam()
summary_stat_ref_closure =  CatalogSummaryStatistics()

function is_valid(param_vector::Vector{Float64})
    global sim_param_closure
    update_sim_param_from_vector!(param_vector,sim_param_closure)
    const rate_tab::Array{Float64,2} = get_any(sim_param_closure, "obs_par", Array{Float64,2})
    #const lambda = sum_kbn(rate_tab)
    if any(x -> x < 0., rate_tab[2:size(rate_tab,1),:]) || any(x -> x > 3., rate_tab[1,:])# || lambda > 10.
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

# https://en.wikipedia.org/wiki/Trigamma_function
function trigamma_x_gr_4(x::T) where T<: Real
   1/x + 0.5/x^2 + 1/(6*x^3) - 1/(30*x^5) + 1/(42*x^7) - 1/(30*x^9) + 5/(66*x^11) - 691/(2730*x^13) + 7/(6*x^15)
end

function trigamma_x_lt_4(x::T) where T<: Real
  n = floor(Int64,5-x)
  z = x+n 
  val = trigamma_x_gr_4(z)
  for i in 1:n
    z -= 1
    val += 1/z^2
  end
  val 
end

function trigamma(x::T) where T<: Real
   x >= 4 ? trigamma_x_gr_4(x) : trigamma_x_lt_4(x)
end


function make_proposal_dist_multidim_beta(theta::AbstractArray{Float64,2}, weights::AbstractArray{Float64,1},  tau_factor::Float64; verbose::Bool = false)
    global sim_param_closure
    const p_dim = length(get_any(sim_param_closure, "p_lim_arr", Array{Float64,1}))-1
    const r_dim = length(get_any(sim_param_closure, "r_lim_arr", Array{Float64,1}))-1
    const max_col_rate = 3.0

    function mom_alpha(x_bar::T, v_bar::T) where T<: Real 
        x_bar * (((x_bar * (1 - x_bar)) / v_bar) - 1)
    end
    function mom_beta(x_bar::T, v_bar::T) where T<: Real 
        (1 - x_bar) * (((x_bar * (1 - x_bar)) / v_bar) - 1)
    end
    # For algorithm, see https://scholarsarchive.byu.edu/cgi/viewcontent.cgi?article=2613&context=etd 
    function fit_beta_mle(x::AbstractArray{T,1}; tol::T = 1e-6, max_it::Int64 = 10, init_guess::AbstractArray{T,1} = Array{T}(undef,0), w::AbstractArray{T,1} = Array{T}(undef,0), verbose::Bool = false ) where T<: Real
        lnxbar =   length(w)>1 ? Compat.Statistics.mean(log.(x),AnalyticWeights(w)) : Compat.Statistics.mean(log.(x))
        ln1mxbar = length(w)>1 ? Compat.Statistics.mean(log.(1.0.-x),AnalyticWeights(w)) : Compat.Statistics.mean(log.(1.0.-x))

        function itterate( mle_guess::Vector{T} ) where T<:Real
            (alpha, beta) = (mle_guess[1], mle_guess[2])
            dgab = digamma(alpha+beta)
            g1 = dgab - digamma(alpha) + lnxbar
            g2 = dgab - digamma(beta) + ln1mxbar
            tgab = trigamma(alpha+beta)
            G = [dgab-trigamma(alpha) tgab; tgab tgab-trigamma(beta)]
            mle_guess -= G \ [g1, g2]
        end 
  
        local mle_new 
        if length(init_guess) != 2
            xbar = length(w)>1 ? Compat.Statistics.mean(x,AnalyticWeights(w)) : Compat.Statistics.mean(x)
            vbar = length(w)>1 ? Compat.Statistics.varm(x,xbar,AnalyticWeights(w)) : Compat.Statistics.varm(x,xbar)
            mle_new = (vbar < xbar*(1.0-xbar)) ? [mom_alpha(xbar, vbar), mom_beta(xbar,vbar)] : ones(T,2)
        else
            mle_new = init_guess
        end
        if verbose
            println("it = 0: ", mle_new)
        end
        if any(mle_new.<=zero(T))
            println("# Warning: mean= ", xbar, " var= ",var," (alpha,beta)_init= ",mle_new," invalid, reinitializing to (1,1)")
            verbose = true
            mle_new = ones(T,2)
        end
        for i in 1:max_it
            mle_old = mle_new
            mle_new = itterate( mle_old )
            epsilon = max(abs.(mle_old.-mle_new))
            if verbose
                println("# it = ", i, ": ", mle_new, " max(Delta alpha, Delta beta)= ", epsilon)
            end
            if epsilon < tol
                break
            end
        end
        return mle_new
    end
    function make_beta(x::AbstractArray{T,1}, w::AbstractArray{T,1}; 
                       mean::T = Compat.Statistics.mean(x,AnalyticWeights(w)), 
                       var::T = Compat.Statistics.varm(x,xbar,AnalyticWeights(w)) ) where T<:Real
        alpha_beta = (var < mean*(1.0-mean)) ? [mom_alpha(mean, var), mom_beta(mean,var)] : ones(T,2)
        if any(alpha_beta.<=zero(T))
            alpha_beta = fit_beta_mle(x, w=w, init_guess=alpha_beta, verbose=true)
        end
        if any(alpha_beta.<=zero(T))
            alpha_beta = ones(T,2)
        end
        Beta(alpha_beta[1], alpha_beta[2])
    end
    function make_beta_transformed(x::AbstractArray{T,1}, w::AbstractArray{T,1}; xmin::T=zero(T), xmax::T=one(T), mean::T = Compat.Statistics.mean(x,AnalyticWeights(w)), var::T = Compat.Statistics.varm(x,xbar,AnalyticWeights(w)) ) where T<:Real
        alpha_beta = (var < mean*(1.0-mean)) ? [mom_alpha(mean, var), mom_beta(mean,var)] : ones(T,2)
        if any(alpha_beta.<=zero(T))
            alpha_beta = fit_beta_mle(x, w=w, init_guess=alpha_beta, verbose=true)
        end
        if any(alpha_beta.<=zero(T))
            alpha_beta = ones(T,2)
        end
        LinearTransformedBeta(alpha_beta[1], alpha_beta[2], xmin=xmin, xmax=xmax)
    end
    
    theta_mean =  sum(theta.*weights',2) # weighted mean for parameters
    tau = tau_factor*ABC.var_weighted(theta'.-theta_mean',weights)  # scaled, weighted covar for parameters
    
    #=
    println("mean= ",theta_mean)
    println("var= ",tau)
    for i in 1:length(theta_mean)
        println("a= ",alpha(theta_mean[i],tau[i]), "  b= ",beta(theta_mean[i],tau[i]))
    end
    =#

    dist_arr = ContinuousDistribution[]
    for j in 1:p_dim
        col_startidx = (j-1)*(r_dim+1)+1
        dist_arr = vcat(dist_arr, make_beta_transformed(theta[col_startidx,:], weights, xmin=0.0, xmax=max_col_rate, mean=theta_mean[col_startidx]/max_col_rate, var=tau[col_startidx]/max_col_rate^2), ContinuousDistribution[ make_beta(theta[i,:], weights, mean=theta_mean[i], var=tau[i]) for i in (col_startidx+1):(col_startidx+r_dim)]   )
    end

    dist = CompositeDist(dist_arr)
end

function make_proposal_dist_multidim_beta(pop::abc_population_type, tau_factor::Float64; verbose::Bool = false)
    make_proposal_dist_multidim_beta(pop.theta, pop.weights, tau_factor, verbose=verbose)
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
    #add_param_fixed(sim_param_closure,"num_kepler_targets",150000)  # For "observed" catalog
    #cat_obs = simulated_read_kepler_observations(sim_param_closure)
    ###
    
    ### Use real planet candidate catalog data
    df_star = setup_star_table_christiansen(sim_param_closure)
    println("# Finished reading in stellar data")
    df_koi,usable_koi = read_koi_catalog(sim_param_closure)
    println("# Finished reading in KOI data")  
    cat_obs = setup_actual_planet_candidate_catalog(df_star, df_koi, usable_koi, sim_param_closure)
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
export setup_abc, run_abc, run_abc_largegen
import ABC
import Distributions
using CompositeDistributions
using Compat
import ExoplanetsSysSim
import EvalSysSimModel
include(joinpath(Pkg.dir(),"ExoplanetsSysSim","examples","poisson_tests", "christiansen_func.jl"))

function setup_abc(num_dist::Integer = 0)
    EvalSysSimModel.setup()
    theta_true = EvalSysSimModel.get_param_vector()
    #param_prior = CompositeDist( Distributions.ContinuousDistribution[Distributions.Uniform(0., 0.3) for x in 1:length(theta_true)] )
    limitP::Array{Float64,1} = get_any(EvalSysSimModel.sim_param_closure, "p_lim_arr", Array{Float64,1})
    const r_dim = length(get_any(EvalSysSimModel.sim_param_closure, "r_lim_arr", Array{Float64,1}))-1
    prior_arr = ContinuousDistribution[]
    for i in 1:(length(limitP)-1)
        max_in_col = floor(3*log(limitP[i+1]/limitP[i])/log(2))
        lambda_col = Uniform(0.0, max_in_col)
        dirch = Dirichlet(ones(r_dim))
        prior_arr = vcat(prior_arr, [lambda_col, dirch])
    end
    param_prior = CompositeDist(prior_arr)
    in_parallel = nworkers() > 1 ? true : false
    
    calc_distance_ltd(sum_stat_obs::ExoplanetsSysSim.CatalogSummaryStatistics,sum_stat_sim::ExoplanetsSysSim.CatalogSummaryStatistics) = EvalSysSimModel.calc_distance(sum_stat_obs,sum_stat_sim,num_dist)
    
    global abc_plan = ABC.abc_pmc_plan_type(EvalSysSimModel.gen_data,EvalSysSimModel.calc_summary_stats,calc_distance_ltd, param_prior, make_proposal_dist=EvalSysSimModel.make_proposal_dist_multidim_beta, is_valid=EvalSysSimModel.is_valid, normalize=EvalSysSimModel.normalize_dirch, num_part=200, num_max_attempt=50, num_max_times=200, epsilon_init=9.9e99, target_epsilon=1.0e-100, in_parallel=in_parallel, adaptive_quantiles = false, epsilon_reduction_factor=0.9, tau_factor=2.0);
end

function run_abc_largegen(pop::ABC.abc_population_type, ss_true::ExoplanetsSysSim.CatalogSummaryStatistics, epshist_targ::Float64, npart::Integer = 1000, num_dist::Integer = 0)
    sim_param_closure = setup_sim_param_christiansen()
    sim_param_closure = set_test_param(sim_param_closure)
    setup_star_table_christiansen(sim_param_closure)
    EvalSysSimModel.set_simparam_ss(sim_param_closure, ss_true)	
    
    theta_true = EvalSysSimModel.get_param_vector()
    #param_prior = CompositeDist( Distributions.ContinuousDistribution[Distributions.Uniform(0., 0.3) for x in 1:length(theta_true)] )
    limitP::Array{Float64,1} = get_any(EvalSysSimModel.sim_param_closure, "p_lim_arr", Array{Float64,1})
    const r_dim = length(get_any(EvalSysSimModel.sim_param_closure, "r_lim_arr", Array{Float64,1}))-1
    prior_arr = ContinuousDistribution[]
    for i in 1:(length(limitP)-1)
        max_in_col = floor(3*log(limitP[i+1]/limitP[i])/log(2))
        lambda_col = Uniform(0.0, max_in_col)
        dirch = Dirichlet(ones(r_dim))
        prior_arr = vcat(prior_arr, [lambda_col, dirch])
    end
    param_prior = CompositeDist(prior_arr)
    in_parallel = nworkers() > 1 ? true : false
    
    calc_distance_ltd(sum_stat_obs::ExoplanetsSysSim.CatalogSummaryStatistics,sum_stat_sim::ExoplanetsSysSim.CatalogSummaryStatistics) = EvalSysSimModel.calc_distance(sum_stat_obs,sum_stat_sim,num_dist)
    
    global abc_plan = ABC.abc_pmc_plan_type(EvalSysSimModel.gen_data,EvalSysSimModel.calc_summary_stats,calc_distance_ltd, param_prior, make_proposal_dist=EvalSysSimModel.make_proposal_dist_multidim_beta, is_valid=EvalSysSimModel.is_valid, normalize=EvalSysSimModel.normalize_dirch, num_part=npart, num_max_attempt=50, num_max_times=1, epsilon_init=9.9e99, target_epsilon=1.0e-100, in_parallel=in_parallel);

    println("# run_abc_largegen: ",EvalSysSimModel.sim_param_closure)
    #if (std(pop.theta, 2)[1])/(mean(pop.theta, 2)[1]) > 0.3
    #    sampler_largegen = abc_plan.make_proposal_dist(pop, 1.2)
    #elseif (std(pop.theta, 2)[1])/(mean(pop.theta, 2)[1]) > 0.1
    #    sampler_largegen = abc_plan.make_proposal_dist(pop, 1.6)
    #else
    sampler_largegen = abc_plan.make_proposal_dist(pop, abc_plan.tau_factor)
    #end
    theta_largegen = Array{Float64}(size(pop.theta, 1), npart)
    weight_largegen = Array{Float64}(npart)
    for i in 1:npart
        theta_val, dist_largegen, attempts_largegen = ABC.generate_theta(abc_plan, sampler_largegen, ss_true, epshist_targ)
        theta_largegen[:,i] = theta_val  
        prior_logpdf = Distributions.logpdf(abc_plan.prior,theta_val)
        sampler_logpdf = ABC.logpdf(sampler_largegen, theta_val)
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
