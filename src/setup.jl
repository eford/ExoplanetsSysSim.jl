## ExoplanetsSysSim/src/setup.jl
## (c) 2015 Eric B. Ford


if VERSION >=v"0.7.0-"
   using Pkg
end

#=
# Is this still needed, now that we've moved to github?
# If so, this would need to be updated for Julia v0.7.0
# Since bitbucket messed up capitalization of package name
if ! isdir( joinpath(Pkg.dir(),"ExoplanetsSysSim") )
     symlink( joinpath(Pkg.dir(),"exoplanetssyssim"), joinpath(Pkg.dir(),"ExoplanetsSysSim") )
end
=#

startdir = pwd()
# Code to be run once to install unregistered package dependencies
# Code to be run once to install non-standard dependencies
try
  Pkg.clone("git@github.com:eford/ABC.jl.git")
  # switch to dragozzine's fork for any minor SysSim specific ABC changes
catch
  warn("Attempted to install ABC.jl package, but was not successful.")
  warn("While most of SysSim will still work, some functionality will not be avaliable unless you install CORBITS correctly.")
end

try
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
catch
  warn("Attempted to install CORBITS.jl package, but was not successful.")
  warn("While most of SysSim will still work, some functionality will not be avaliable unless you install CORBITS correctly.")
end

