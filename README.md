# ExoplanetsSysSim
The package does not yet work.
Currently, we're implementing a generic framework.

# How to get started:
* Create bitbucket.org account, get permission to access this repository
* Make sure you have julia and git installed, and an ssh key
* Person Icon (up, right), Manage Account, SSH keys, Add Key.  Paste contents of cat ~/.ssh/id_rsa.pub into box. Ok. Also need to set this up for GitHub because you'll use Eric's ABC package there.
* Install this repo as a package and it's non-standard dependencies (ABC, CORBITS)
```
#!julia

   julia
   Pkg.clone("git@bitbucket.org:eford/exoplanetssyssim.git")
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
~/.julia/v0.4/ExoplanetsSysSim (or the similar appropriate path).

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

