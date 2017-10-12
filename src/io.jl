## ExoplanetsSysSim/src/io.jl
## (c) 2015 Eric B. Ford

module SysSimIO
using ExoplanetsSysSim
using HDF5, JLD

#if VERSION >= v"0.5-"
#  import Compat: UTF8String, ASCIIString
#end

export save_sim_param, save_sim_results, load_sim_param, load_distances, load_summary_stats

function save_sim_param(filename::String, p::SimParam)
 local file
 try
   file = JLD.jldopen(filename,"w")
   write_sim_param(file,p)
 catch
   println("# Problem writing parameters to jld file: ", filename)
   return false
 finally
   close(file)
 end
 return true
end

function save_sim_results(filename::String, p::SimParam; distances::Vector{Float64}=Array{Float64}(0), summary_stats::CatalogSummaryStatistics = CatalogSummaryStatistics() )
 local file
 try
   file = JLD.jldopen(filename,"w")
   write_sim_param(file,p)
   write_sim_summary_stats(file,summary_stats)
   write_sim_distances(file,distances)
 catch
   println("# Problem writing data to jld file: ", filename)
 finally
   close(file)
  end
end

function write_sim_distances(file::JLD.JldFile, d::Vector{Float64} )
  JLD.write(file,"distances",d)
end

function write_sim_summary_stats(file::JLD.JldFile, ss::CatalogSummaryStatistics )
  JLD.write(file,"summary_stats",ss.stat)
end

function write_sim_param(file::JLD.JldFile, p::SimParam)
 sim_param_bool = Dict{String,Bool}()
 sim_param_int = Dict{String,Integer}()
 sim_param_real = Dict{String,Real}()
 sim_param_function = Dict{String,String}()
 sim_param_string = Dict{String,String}()
 for k in keys(p.param)
   #println("# k=",k,".  v=",p.param[k])
   if typeof(p.param[k]) <: Bool
      sim_param_bool[k] = p.param[k]
   elseif typeof(p.param[k]) <: Integer
      sim_param_int[k] = p.param[k]
   elseif typeof(p.param[k]) <: Real
      sim_param_real[k] = p.param[k]
   elseif typeof(p.param[k]) <: Function
      sim_param_function[k] = string(p.param[k])
   elseif typeof(p.param[k]) <: AbstractString
      sim_param_string[k] = convert(String,p.param[k])
   else
	  warn(string("Can't store value of >",k,"< due to type ", typeof(p.param[k])))
   end
 end
   JLD.write(file,"sim_param_int",sim_param_int)
   JLD.write(file,"sim_param_real",sim_param_real)
   JLD.write(file,"sim_param_function",sim_param_function)
   JLD.write(file,"sim_param_string",sim_param_string)
   JLD.write(file,"sim_param_bool",sim_param_bool)
   JLD.write(file,"sim_param_active",p.active)
end


function load_sim_param(filename::String)
  local jld_data
  try  
    jld_data = load(filename)
  catch
   println("# Problem reading parameters from jld file: ", filename) 
  end
  p = SimParam()
  merge!(p.active, jld_data["sim_param_active"])
  merge!(p.param, jld_data["sim_param_int"])
  merge!(p.param, jld_data["sim_param_string"])
  merge!(p.param, jld_data["sim_param_real"])
  #=
   merge!(p.param, jld_data["sim_param_function"]) 
  df = jld_data["sim_param_function"]
  for k in keys(df)
    p.param[k]::Function = symbol(df[k])
  end
  =#
  merge!(p.param, jld_data["sim_param_bool"])
  return p
end

function load_distances(filename::String)
  local jld_data
  try  
    jld_data = load(filename)
  catch
   println("# Problem reading distances from jld file: ", filename) 
  end
  d::Array{Float64,1} = jld_data["distances"]
  return d
end

function load_summary_stats(filename::String)
  local jld_data
  try  
    jld_data = load(filename)
  catch
   println("# Problem reading parameters from jld file: ", filename) 
  end
  s = CatalogSummaryStatistics()
  merge!(s.stat, jld_data["summary_stats"])
  return s
end

function test_io()
  sim_param = setup_sim_param_demo()
  save_sim_param("test.jld",sim_param)
  spd = load_sim_param("test.jld")
  rm("test.jld")
end

end # module SysSimIO
