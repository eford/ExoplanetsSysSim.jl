# SysSim 

# To-Do List, short-term
  * Add more flexible period-size distribution to planetary_system.jl
  * Generate systems with orbital planes clustered around fundamental plane
  * Add summary statistics that we know we'll want to summary_stats.jl
  * Add distances that we know we'll want to abc_distances.jl
  * Add more tests

Determine which (if any) of the following are important to add for the first paper:
  * measurements of host star properties (What's a good enough approximation?  We should test if our results will be sensitive to this, as we might be able to get away with something simplistic like truncated Normals), 
  * multiple CDPPs for different quarters/seasons and transit durations (wouldn't be hard, but substantially adds to memory consumption),
  * contamination (can we just look up values from some other catalog?),
  * interpolate MES threshold for different transit durations (I question whether this is significant),
  * planets/BDs/EBs around secondary/background stars (ask Robert for help?), 
  * limb darkening (currently have an rough approximation for a correction factor in, so probalby good enough), 
  * other complexities for transit detection probability model (e.g., actual data gaps, short cadence)? 
  * [add more here],
  *  etc.

# To-Do List, medium-term
  * Make sure algorithm can perform well on simulated data
  * Perform initial analysis of existing catalog
  * Write first paper

# Before using for any Science
  * grep IMPORTANT *.jl and solve or demote those issues
  * grep WARNING *.jl and be aware of these limitationrs (or resolve them)
  * grep TODO *.jl and be aware of these limitationrs (or resolve them)
  * grep TODO *.jl |grep -v OPT | grep -v DETAIL to skip the ones that probably aren't critical for the first paper or are just suggestions on how we could make the code more efficient
  * grep QUERY *.jl for long-term issues and discuss what to do about them.


# To-Do List, long-term

# Notes to develoeprs:
  * Some files that just place holders (i.e., not yet used): koi_table.jl, limb_darkening.jl


