## ExoplanetsSysSim/examples/run_abc_parallel.jl
## (c) 2015 Eric B. Ford

if nworkers() == 1
  addprocs(8)
end

import ExoplanetsSysSim
import ABC
@everywhere include("abc_setup.jl")

import SysSimABC; @everywhere using SysSimABC

@everywhere abc_plan = setup_abc(3)

@time output = SysSimABC.run_abc(abc_plan)

@time output2 = SysSimABC.run_abc(abc_plan,output)
#@time output3 = SysSimABC.run_abc(abc_plan,output2)

