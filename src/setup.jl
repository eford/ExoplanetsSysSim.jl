## ExoplanetsSysSim/src/setup.jl
## (c) 2015 Eric B. Ford

# Since bitbucket messed up capitalization of package name
if ! isdir( joinpath(Pkg.dir(),"ExoplanetsSysSim") )
     symlink( joinpath(Pkg.dir(),"exoplanetssyssim"), joinpath(Pkg.dir(),"ExoplanetsSysSim") )
end

# Code to be run once to install unregistered package dependencies
# Code to be run once to install non-standard dependencies
Pkg.clone("git@github.com:eford/ABC.jl.git")
# switch to dragozzine's fork for any minor SysSim specific ABC changes

Pkg.clone("git@github.com:jbrakensiek/CORBITS.git")

# Compile CORBITS library and put it somewhere we can find
cd(joinpath(Pkg.dir(),"CORBITS"))
run(`make lib`)
cd(homedir())
if !is_windows()
   symlink( joinpath(Pkg.dir("CORBITS"),"libcorbits.so"), joinpath(Pkg.dir("ExoplanetsSysSim"),"libcorbits.so") )
else
   cp( joinpath(Pkg.dir("CORBITS"),"libcorbits.so"), joinpath(Pkg.dir("ExoplanetsSysSim"),"libcorbits.so") )
end


