using JLD
using ExoplanetsSysSim
using ABC
using StatsBase

dir = "hz320_q1q16"
limitP = [237., 320.]
limitR = [1.0, 1.5]
dens_denom = 1.0/(log2(limitP[2])-log2(limitP[1]))/(log2(limitR[2])-log2(limitR[1]))

file_out = open("mc_rate_table.txt", "w")

println("Reading in ", dir)
pop_out = load(string("./",dir,"/test-1-out.jld"), "output")

weight_vec = pweights(pop_out.weights)
quant_arr = quantile(pop_out.theta[1,:], weight_vec, [0.1587, 0.5, 0.8413])
write(file_out, string(quant_arr[2], " + ", quant_arr[3]-quant_arr[2], " - ", quant_arr[2]-quant_arr[1], "\n"))
write(file_out, string(quant_arr[2]*dens_denom, " + ", (quant_arr[3]-quant_arr[2])*dens_denom, " - ", (quant_arr[2]-quant_arr[1])*dens_denom, "\n"))

close(file_out)
