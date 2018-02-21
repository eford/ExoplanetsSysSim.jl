# ExoplanetsSysSim
The package is still considered experimental, while we work to implement a generic framework.

# How to get started:
* Make sure you have julia and git installed.
* If you want to use ssh keys instead of https authentication, then:
  * setup a local ssh key
  * Tell Github about your ssh key:  Person Icon (upper right), Settings, SSH & GPG keys, New SSH Key.  Entry a name in the title box and paste the contents of `cat ~/.ssh/id_rsa.pub` into the "Key" box. Add SSH Key.  
* Install the ExoplanetsSysSim repo as a Julia package using the following julia command.  The setup script below will attempt to install two related packages (ABC, CORBITS) that are often used in combination with ExoplanetsSysSim.
```
#!julia

   #julia
   Pkg.clone("git@github.com:eford/ExoplanetsSysSim.jl.git")
   include(joinpath(Pkg.dir(),"ExoplanetsSysSim/src/setup.jl"))   

```
* If you are using windows, you might encoutner an issue with capitalization of package names.
* If you have some issues with Blosc, you might need to follow the instructions here: https://github.com/stevengj/Blosc.jl/issues/5
* Run some tests, e.g. 
```
#!julia

   using ExoplanetsSysSim
   include(joinpath(Pkg.dir(),"ExoplanetsSysSim/examples/basic/test.jl"))   
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
  * Matthias He
  * Danley Hsu
  * Darin Ragozzine
# Other Contributors/Consultants:
  * Robert Morehead
  * Keir Ashby
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
