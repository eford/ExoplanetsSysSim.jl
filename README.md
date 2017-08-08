# ExoplanetsSysSim
The package is still considered experimental, while we work to implement a generic framework.

# How to get started:
* Create github.com account, and (until it's made public) get permission to access this repository
* Make sure you have julia and git installed, and an ssh key
* Person Icon (upper right), Settings, SSH & GPG keys, New SSH Key.  Entry a name in the title box and paste the contents of `cat ~/.ssh/id_rsa.pub` into the "Key" box. Add SSH Key.  
* Install this repo as a package and it's non-standard dependencies (ABC, CORBITS)
```
#!julia

   julia
   Pkg.clone("git@github.com:eford/ExoplanetsSysSim.jl.git")
   include(joinpath(Pkg.dir(),"exoplanetssyssim/src/setup.jl"))   

```
* If you have some issues with Blosc, you might need to follow the instructions here: https://github.com/stevengj/Blosc.jl/issues/5
* Run some tests, e.g. 
```
#!julia

   using ExoplanetsSysSim
   include(joinpath(Pkg.dir(),"exoplanetssyssim/examples/test.jl"))   
```
* Run some simple "applications" (after exiting out of Julia and going to 
~/.julia/v0.6/ExoplanetsSysSim (or the similar appropriate path).

```
#!csh
cd apps
julia syssim_summary_stats.jl demo_param.in demo_ss.out
julia syssim_dist.jl demo_param.in demo_ss.out
```
* Create your own feature branch and start making SysSim more realistic

# Team:
# Developers:
  * Eric Ford
  * Danley Hsu
  * Darin Ragozzine
# Other Contributors/Consultants:
  * Robert Morehead
  * Jessi Cisewski
  * Chad Schafer
  * Tom Loredo
  * Robert Wolpert

# Acknowledgements:
* NASA
  * Kepler Mission
  * Kepler Science Team
  * Kepler Multi-body & Transit Timing Variations Working Groups
  * Origins of Solar Systems program, award NNX14AI76G
* The Pennsylvanis State University
  * Dept. of Astronomy & Astrophysics
  * Center for Exoplanets & Habitable Worlds
  * Eberly College of Science
  * Institute for CyberScience
  * Center for Astrostatistics
  * Penn State Astrobiology Research Center
* Florida Institute of Technology
* University of Florida
* Statistical and Applied Mathematical Sciences Institute
