## ExoplanetsSysSim/examples/run_abc.jl
## (c) 2015 Eric B. Ford

include("abc_setup.jl")

using SysSimABC

abc_plan = setup_abc(3)
@time output = SysSimABC.run_abc(abc_plan)

#@time output2 = SysSimABC.run_abc(abc_plan,output)
#@time output3 = SysSimABC.run_abc(abc_plan,output2)

