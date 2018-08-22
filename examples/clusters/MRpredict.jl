if !isdefined(:ExoplanetsSysSim) using ExoplanetsSysSim end
using Roots #if we want to use the 'find_zero()' function
#using Optim

##### This code is a translation of 'MRpredict.R' at 'https://github.com/Bo-Ning/Predicting-exoplanet-mass-and-radius-relationship/blob/master/MR-predict/MRpredict.R'

struct MR_param_Ning2018
    Mass_min::Float64 #log10 scale
    Mass_max::Float64 #log10 scale
    Radius_min::Float64 #log10 scale
    Radius_max::Float64 #log10 scale
    deg::Real
    weights_mle::Array{Float64,1}
end

function pdfnorm_beta(x, x_obs, x_sd, x_max::Float64, x_min::Float64, shape1::Real, shape2::Real ; log::Bool=true)
    if log == true
        norm_dist = Normal(10^x, x_sd)
    else
        norm_dist = Normal(x, x_sd)
    end
    beta_dist = Beta(shape1, shape2)
    norm_beta = pdf.(norm_dist, x_obs) .* pdf.(beta_dist, (x-x_min)/(x_max-x_min)) / (x_max-x_min)
    return norm_beta
end

function fn_for_integrate(data, deg::Real, degree::Real, x_max::Float64, x_min::Float64 ; log::Bool=false, abs_tol::Float64=1e-10)
    x_obs, x_sd = data[1], data[2]
    shape1, shape2 = degree, deg-degree+1
    return quadgk(x -> pdfnorm_beta(x, x_obs, x_sd, x_max, x_min, shape1, shape2; log=log), x_min, x_max, abstol=abs_tol)[1]
end

function conditional_density(y::Real, y_max::Float64, y_min::Float64, x_max::Float64, x_min::Float64, deg::Real, w_hat::Vector{Float64} ; y_sd=nothing, qtl::Vector{Float64}=[0.16, 0.84]) ###is 'y_sd=nothing' a good choice?
    deg_vec = 1:deg

    #To compute conditional mean, variance, quantile, distribution:
    if y_sd == nothing
        y_stdardize = (y-y_min)/(y_max-y_min)
        y_beta_indv = map(degree -> pdf.(Beta(degree, deg-degree+1), y_stdardize), deg_vec) / (y_max-y_min)
    else
        y_beta_indv = map(degree -> fn_for_integrate([y, y_sd], deg, degree, y_max, y_min), deg_vec)
    end
    y_beta_pdf = kron(ones(deg), y_beta_indv)
    denominator = sum(w_hat .* y_beta_pdf)

    #Mean:
    mean_beta_indv = deg_vec/(deg+1)*(x_max-x_min)+x_min
    mean_beta = kron(mean_beta_indv, y_beta_indv)
    mean_numerator = sum(w_hat .* mean_beta)
    mean = mean_numerator/denominator

    #Variance:
    var_beta_indv = deg_vec .* (deg-deg_vec+1) / ((deg+1)^2*(deg+2))*(x_max-x_min)^2
    var_beta = kron(var_beta_indv, y_beta_indv)
    var_numerator = sum(w_hat .* var_beta)
    var = var_numerator/denominator

    #Quantile:
    function pbeta_conditional_density(x)
        function mix_density(j)
            x_indv_cdf = map(degree -> cdf.(Beta(degree, deg-degree+1), (j-x_min)/(x_max-x_min)), deg_vec)
            quantile_numerator = sum(w_hat .* kron(x_indv_cdf, y_beta_indv))
            return p_beta = quantile_numerator/denominator
        end

        return map(mix_density, x)
    end

    function conditional_quantile(q)
        function g(x)
            return pbeta_conditional_density(x) - q
        end
        return find_zero(g, (x_min, x_max)) ###'find_zero' is from 'using Roots' (is this a good function to use?); Eric suggested brent's method in the Optim package
        #return Optim.minimum(optimize(g, x_min, x_max))
    end

    quantile = map(conditional_quantile, qtl)

    return [mean, var, quantile, denominator, y_beta_indv]
end

function predict_mass_given_radius(Radius::Real, param::MR_param_Ning2018 ; R_sigma=nothing, posterior_sample::Bool=false, qtl::Vector{Float64}=[0.16, 0.84])
    #Convert data to log scale:
    l_radius = log10.(Radius)

    #The 'posterior_sample == false' condition can deal with two cases:
    #Case I: if input data do not have measurement errors
    #Case II: if input data have measurement error
    if posterior_sample == false
        predicted_value = conditional_density(l_radius, param.Radius_max, param.Radius_min, param.Mass_max, param.Mass_min, param.deg, param.weights_mle; y_sd=R_sigma, qtl=qtl)
        predicted_mean = predicted_value[1] ###
        predicted_lower_quantile = predicted_value[3][1] ###
        predicted_upper_quantile = predicted_value[3][2] ###
    elseif posterior_sample == true
        #Case III: if the input are posterior samples
        radius_sample = log10(Radius)
        k = length(radius_sample)
        denominator_sample = zeros(k)
        mean_sample = zeros(k)
        y_beta_indv_sample = zeros(k, param.deg)

        #Calculate mean:
        for i in 1:k
            results = cond_density_estimation(y=radius_sample[i], y_max=param.Radius_max, y_min=param.Radius_min, x_max=param.Mass_max, x_min=param.Mass_min, deg=param.deg, w_hat=param.weights_mle, qtl=quantile, only_output_mean=true) ###what function is this??? Is it part of some CRAN package...
            mean_sample[i] = (results[1]) ###why are there parentheses here? It makes no difference...
            denominator_sample[i] = results[2]
            y_beta_indv_sample[i,:] = results[3:57]
        end
        predicted_mean = mean(mean_sample)

        #Calculate 16% and 84% quantiles:
        #Mixture of the CDF of k conditional densities
        function pbeta_conditional_density(x)
            function mix_density(j)
                deg_vec = 1:param.deg
                x_indv_cdf = map(degree -> cdf.(Beta(degree, param.deg-degree+1), (j-x_min)/(x_max-x_min)), deg_vec)
                quantile_numerator = zeros(k)
                p_beta_sample = zeros(k)
                for ii in 1:k
                    quantile_numerator[ii] = sum(param.weights_mle .* kron(x_indv_cdf, y_beta_indv_sample[ii,:]))
                    p_beta_sample[ii] = quantile_numerator[ii]/denominator_sample[ii]
                end
                return p_beta = mean(p_beta_sample)
            end

            return map(mix_density, x)
        end

        function mixture_conditional_quantile(q, x_min, x_max)
            function g(x)
                return pbeta_conditional_density(x) - q
            end
            return find_zero(g, (x_min, x_max)) ###see note above
            #return Optim.minimum(optimize(g, x_min, x_max))
        end

        predicted_quantiles = map(q -> mixture_conditional_quantile(q, param.Mass_min, param.Mass_max), qtl)
        predicted_lower_quantile = predicted_quantiles[1]
        predicted_upper_quantile = predicted_quantiles[2]
    end

    #Return the output:
    return [predicted_mean, predicted_lower_quantile, predicted_upper_quantile]
end

function draw_planet_mass_from_radius_Ning2018(Radius::Real, param::MR_param_Ning2018)
    #This function takes in a Radius (in solar radii) and draws a mass (returning in solar masses) probabilistically
    @assert param.Radius_min < log10.(Radius/ExoplanetsSysSim.earth_radius) < param.Radius_max
    q = rand()
    l_mass = predict_mass_given_radius(Radius, param; qtl=[q,q])[2]
    return mass = (10^l_mass)*ExoplanetsSysSim.earth_mass
end

function generate_planet_mass_from_radius_Ning2018(Radius::Real, sim_param::SimParam)
    global MR_param
    return draw_planet_mass_from_radius_Ning2018(Radius, MR_param)
end





##### Examples:

weights_mle = CSV.read("weights.mle.csv")[:x]
MR_param = MR_param_Ning2018(-1., 3.809597, -0.302, 1.357509, 55, weights_mle) #(-1., 3.809597, -0.3, 1.357509, 55, weights_mle)

#Observation without measurement errors:
#=
Radius = 5 #original scale, not log scale
predict_result = predict_mass_given_radius(Radius, MR_param)
println(predict_result)
=#

#Observation with a measurement error:
#=
Radius = 5 #original scale, not log scale
R_sigma = 0.1
predict_result = predict_mass_given_radius(Radius, MR_param; R_sigma=0.1)
println(predict_result)
=#

#Input are posterior samples: ###currently broken because the function 'cond_density_estimation' is undefined
#=
Radius_dist = Normal(5, 0.5) #original scale, not log scale
Radius_samples = rand(Radius_dist, 100)
predict_result = predict_mass_given_radius(Radius_samples, MR_param; R_sigma=0.1, posterior_sample=true)
println(predict_result)
=#

#If want to change 16% and 84% quantiles to 5% and 95% quantiles:
#=
Radius = 5 #original scale, not log scale
predict_result = predict_mass_given_radius(Radius, MR_param; qtl=[0.05, 0.95])
println(predict_result)
=#

