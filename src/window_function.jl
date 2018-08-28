## ExoplanetsSysSim/src/window_function.jl
## (c) 2018 Darin Ragozzine

# Gather and prepare the window function data
module WindowFunction

export setup_window_function, get_window_function_data, get_window_function_id, eval_window_function

#using DataArrays
using DataFrames
#using CSV
using JLD2
using ExoplanetsSysSim.SimulationParameters


# Object to hold window function data
struct window_function_data
  window_func_array::Array{Float64,3}      # Value of window function (window_function_id, duration_id, period_id).  Maybe rename to wf_value or data?
  wf_durations_in_hrs::Array{Float64,1}    # Boundaries for duration bins in window_func_array.  Maybe rename to durations?
  wf_periods_in_days::Array{Float64,1}     # Boundaries for periods bins in window_func_array.  Maybe rename to periods?
  sorted_quarter_strings::Array{Int64,1}   # TODO OPT: Is there a reason to keep this?  Maybe rename to quarter_strings?
  allsortedkepids::Array{Int64,1}          # value is Kepler ID.  Index is same as index to window_function_id_arr
  window_function_id_arr::Array{Int64,1}   # value is index for window_func_array.  Index is same as index to allsortedkepids
  default_wf_id::Int64                     # Index corresponding to the default window function
end


function window_function_data()
  window_function_data( Array{Float64,3}(undef,0,0,0), Array{Float64,1}(undef,0),Array{Float64,1}(undef,0), Array{Int64,1}(undef,0),Array{Int64,1}(undef,0),Array{Int64,1}(undef,0), 0 )
end


win_func_data = window_function_data()

function setup(sim_param::SimParam; force_reread::Bool = false)
  global win_func_data
  if haskey(sim_param,"read_window_function") && !force_reread
     return win_func_data
  end
  window_function_filename = convert(String,joinpath(dirname(pathof(ExoplanetsSysSim)), "data", convert(String,get(sim_param,"window_function","DR25topwinfuncs.jld2")) ) )
  setup(window_function_filename)
  add_param_fixed(sim_param,"read_window_function",true)
  @assert( size(win_func_data.window_func_array,2) == length(win_func_data.wf_durations_in_hrs) )
  @assert( size(win_func_data.window_func_array,3) == length(win_func_data.wf_periods_in_days) )
  @assert( size(win_func_data.window_func_array,1) >= maximum(win_func_data.window_function_id_arr) )
  @assert( size(win_func_data.window_func_array,1) >= win_func_data.default_wf_id )
  return win_func_data
end



function setup(filename::String)
# Reads in the window function data collected from the Kepler Completeness Products
# see Darin Ragozzine's get/cleanDR25winfuncs.jl

  if ismatch(r".jld2$",filename)
    try
      wfdata = load(filename)
      window_func_array = wfdata["window_func_array"]
      wf_durations_in_hrs = wfdata["wf_durations_in_hrs"]  # TODO OPT DETAIL: Should we convert units to days here?
      wf_periods_in_days = wfdata["wf_periods_in_days"]
      sorted_quarter_strings = wfdata["sorted_quarter_strings"]
      allsortedkepids = wfdata["allsortedkepids"]
      window_function_id_arr = wfdata["window_function_id_arr"]

      global win_func_data = window_function_data(window_func_array, wf_durations_in_hrs, wf_periods_in_days, sorted_quarter_strings, 
                                                allsortedkepids, window_function_id_arr, maximum(window_function_id_arr) )

    catch
      error(string("# Failed to read window function data > ", filename," < in jld2 format."))
    end
  end
 
  return win_func_data
end

setup_window_function(sim_param::SimParam; force_reread::Bool = false) = setup(sim_param, force_reread=force_reread)
setup_window_function(filename::String; force_reread::Bool = false) = setup(filename, force_reread=force_reread)

function get_window_function_data()::window_function_data
   #global win_func_data
   return win_func_data
end

function get_window_function_id(kepid::Int64; use_default_for_unknown::Bool = true)::Int64
  # takes the quarter string from the stellar catalog and determines the window function id
  # from DR25topwinfuncs.jld2 made by Darin Ragozzine's cleanDR25winfuncs.jl script.
  no_win_func_available::Int64 = -1        # hardcoding this in, should match convention in window function input file

  wf_id = win_func_data.window_function_id_arr[searchsortedfirst(win_func_data.allsortedkepids,kepid)] # all Kepler kepids are in allsortedkepids

  if wf_id == no_win_func_available && use_default_for_unknown
    # if a target is observed for less than 4 quarters, then it won't have a corresponding
    # window function in this list, so throw a warning and use the last window_function_id
    # which corresponds to an "averaged" window function
    warn("Window function data is not avaialble for kepid $kepid, using default.")
    wf_id = win_func_data.default_wf_id
  end
  # TODO SCI DETAIL IMPORTANT? This does not include TPS timeouts or MESthresholds (see DR25 Completeness Products)

  return wf_id 
end


function calc_period_idx(P::Float64)::Int64
  @assert(P>zero(P))
  idx = searchsortedlast(win_func_data.wf_periods_in_days,P)
  if idx == 0
     return 1
  elseif idx<length(win_func_data.wf_periods_in_days)
     if P-win_func_data.wf_periods_in_days[idx]>win_func_data.wf_periods_in_days[idx+1]-P
        idx += 1
     end 
  end
  return idx   # TODO IMPORTANT: IMPLEMENT / TEST
end

function calc_duration_idx(D::Float64)::Int64 
  # NOTE IMPORTANT: Currently assumes we left wf data in hours, so deal with that conversion here
  @assert(D>zero(D))
  hours_in_day = 24 
  idx = searchsortedlast(win_func_data.wf_durations_in_hrs,D*hours_in_day)
  if idx == 0
     return 1
  elseif idx<length(win_func_data.wf_durations_in_hrs)
     if D*hours_in_day-win_func_data.wf_durations_in_hrs[idx]>win_func_data.wf_durations_in_hrs[idx+1]-D*hours_in_day
        idx += 1
     end
  end
  return idx   # TODO IMPORTANT: IMPLEMENT / TEST
end


function eval_window_function(wf_idx::Int64=-1; Duration::Float64=0., Period::Float64=0.)::Float64
  D_idx = calc_duration_idx(Duration)
  P_idx = calc_period_idx(Period)
  wf = eval_window_function(wf_idx,D_idx,P_idx)
  # TODO IMPORTANT: Improve way deal with missing wf values for some durations. Interpolate?
  while wf<=zero(wf) && D_idx<length(win_func_data.wf_durations_in_hrs)  
     D_idx += 1
     wf = eval_window_function(wf_idx,D_idx,P_idx)
  end
  return wf
end

function eval_window_function(wf_idx::Int64, D_idx::Int64, P_idx::Int64)::Float64
   global win_func_data
   #@assert(1<=wf_idx<maximum(win_func_data.window_function_id_arr))
   #@assert(1<=P_idx<=length(win_func_data.wf_periods_in_days))
   #@assert(1<=D_idx<=length(win_func_data.wf_durations_in_hrs))
   return win_func_data.window_func_array[wf_idx,D_idx,P_idx]     
end


end  # module WindowFunction



