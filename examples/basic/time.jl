## ExoplanetsSysSim/examples/test.jl
## (c) 2015 Eric B. Ford

workspace()
using ExoplanetsSysSim

println("# Starting end-to-end test (setting parameters, creating physical catalogs, creating observed catalogs, computing summary statistics, computing the distance between them).  ")
@time test_eval_model()
@time test_eval_model()

