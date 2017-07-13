## ExoplanetsSysSim/examples/run_abc.jl
## (c) 2015 Eric B. Ford

include(joinpath(Pkg.dir(),"ExoplanetsSysSim","examples","multiplanet_systems", "model.jl"))
include("abc_setup.jl")

using SysSimABC
using JLD

start_ind = 1
final_ind = 1

for n in start_ind:final_ind
  abc_plan = setup_abc()

  param_names = make_vector_of_active_param_keys(EvalSysSimModel.sim_param_closure)
  println("# Param Names = ", param_names)
  @time output = SysSimABC.run_abc(abc_plan)

  println("")
  println("Test ", n)
  println("")
  println("Mean = ")
  println(reshape(mean(output.theta,2), (1, size(output.theta,1))))
  println("")
  println("Standard deviation = ")
  println(reshape(std(output.theta,2), (1, size(output.theta,1))))
  println("")
  #println(EvalSysSimModel.get_ss_obs())
  println("-----------------------------")
  println("# Param Names = ", param_names)
  save(string("test-",n,"-out.jld"), "output", output, "ss_true", EvalSysSimModel.get_ss_obs())
  sample = SysSimABC.run_abc_largegen(output,EvalSysSimModel.get_ss_obs(),output.accept_log.epsilon[end])
  save(string("test-",n,"-sample.jld"), "output", output, "ss_true", EvalSysSimModel.get_ss_obs(), "sample", sample)
end

