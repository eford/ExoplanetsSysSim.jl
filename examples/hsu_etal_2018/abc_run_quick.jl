## ExoplanetsSysSim/examples/abc_run.jl
## (c) 2015 Eric B. Ford

include(joinpath(pwd(), "abc_setup.jl"))

using SysSimABC
using JLD
using StatsBase

limitP = [237., 320.]
limitR = [1.0, 1.5]
dens_denom = 1.0/(log2(limitP[2])-log2(limitP[1]))/(log2(limitR[2])-log2(limitR[1]))

abc_plan = setup_abc(1)
EvalSysSimModel.add_param_fixed(EvalSysSimModel.sim_param_closure,"num_targets_sim_pass_one",1000)  # Set universal simulated catalog size
@time output = SysSimABC.run_abc(abc_plan)

weight_vec = pweights(output.weights)
quant_arr = quantile(output.theta[1,:], weight_vec, [0.1587, 0.5, 0.8413])

println("")
println("-----------------------------")
println("Rate = ", string(quant_arr[2], " + ", quant_arr[3]-quant_arr[2], " - ", quant_arr[2]-quant_arr[1]))
println("")
println("Density = ", string(quant_arr[2]*dens_denom, " + ", (quant_arr[3]-quant_arr[2])*dens_denom, " - ", (quant_arr[2]-quant_arr[1])*dens_denom))
println("")
println(EvalSysSimModel.get_ss_obs())

file_out = open("rate_output.txt", "w")
write(file_out, string(quant_arr[2], " + ", quant_arr[3]-quant_arr[2], " - ", quant_arr[2]-quant_arr[1], "\n"))
write(file_out, string(quant_arr[2]*dens_denom, " + ", (quant_arr[3]-quant_arr[2])*dens_denom, " - ", (quant_arr[2]-quant_arr[1])*dens_denom, "\n"))
close(file_out)

save(string("test-pop-out.jld"), "output", output, "ss_true", EvalSysSimModel.get_ss_obs())
