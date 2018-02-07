#using ExoplanetsSysSim
#using StatsFuns
if !isdefined(:JLD) using JLD end
if !isdefined(:DataFrames) using DataFrames end
import DataFrames.DataFrame, DataFrames.isna
#import ExoplanetsSysSim.StellarTable.df
#import ExoplanetsSysSim.StellarTable.usable
import Compat: UTF8String, ASCIIString

## Old code for generating stellar properties # TODO: WARNING: Should eventually use version in main branch to make sure have recent improvements

## stellar_table
function setup_star_table_christiansen(sim_param::SimParam; force_reread::Bool = false)
  #global df
  df = ExoplanetsSysSim.StellarTable.df
  if haskey(sim_param,"read_stellar_catalog") && !force_reread
     return df
     #return data
  end
  stellar_catalog_filename = convert(ASCIIString,joinpath(Pkg.dir("ExoplanetsSysSim"), "data", convert(ASCIIString,get(sim_param,"stellar_catalog","q1_q17_dr25_stellar.csv")) ) )
  df = setup_star_table_christiansen(stellar_catalog_filename)
  add_param_fixed(sim_param,"read_stellar_catalog",true)
  set_star_table(df)
  return df  
end

function setup_star_table_christiansen(filename::ASCIIString; force_reread::Bool = false)
  #global df, usable
  df = ExoplanetsSysSim.StellarTable.df
  usable = ExoplanetsSysSim.StellarTable.usable
  if ismatch(r".jld$",filename)
  try 
    data = DataFrames.load(filename)
    df::DataFrame = data["stellar_catalog"]
    usable::Array{Int64,1} = data["stellar_catalog_usable"]
    set_star_table(df, usable)
  catch
    error(string("# Failed to read stellar catalog >",filename,"< in jld format."))
  end
  println("Test 1?")
  else
  try 
    df = DataFrames.readtable(filename)
  catch
    error(string("# Failed to read stellar catalog >",filename,"< in ascii format."))
  end
  end # if ismatch
  println("Test 2?")
  has_mass = ! (isna.(df[:mass]) | isna.(df[:mass_err1]) | isna.(df[:mass_err2]))
  has_radius = ! (isna.(df[:radius]) | isna.(df[:radius_err1]) | isna.(df[:radius_err2]))
  has_dens = ! (isna.(df[:dens]) | isna.(df[:dens_err1]) | isna.(df[:dens_err2]))
  has_rest = ! (isna.(df[:rrmscdpp04p5]) | isna.(df[:dataspan]) | isna.(df[:dutycycle]))
  in_Q1Q12 = []
#=for x in df[:st_quarters]
    subx = string(x)
    subx = ("0"^(17-length(subx)))*subx
    indQ = search(subx, '1')
    if ((indQ < 1) | (indQ > 12))
      push!(in_Q1Q12, false)
    else
      push!(in_Q1Q12, true)
    end
  end
=#
  is_FGK = []
  for x in 1:length(df[:teff])
    if ((df[x,:teff] > 4000.0) & (df[x,:teff] < 7000.0) & (df[x,:logg] > 4.0))
#println("logg?: ", df[x,:logg])
      push!(is_FGK, true)
    else
      push!(is_FGK, false)
    end
  end   
  is_usable = has_radius & is_FGK & has_mass & has_rest #& in_Q1Q12 # & has_dens
  #if contains(filename,"q1_q12_christiansen.jld")
#if contains(filename,"q1_q12_christiansen")   # TODO: Ask Danely what he's trying to do here.
#is_usable = is_usable #& in_Q1Q12
#end
  # See options at: http://exoplanetarchive.ipac.caltech.edu/docs/API_keplerstellar_columns.html
  # TODO SCI DETAIL or IMPORTANT?: Read in all CDPP's, so can interpolate?
  symbols_to_keep = [ :kepid, :mass, :mass_err1, :mass_err2, :radius, :radius_err1, :radius_err2, :dens, :dens_err1, :dens_err2, :rrmscdpp04p5, :dataspan, :dutycycle ]
  delete!(df, [~(x in symbols_to_keep) for x in names(df)])    # delete columns that we won't be using anyway
  usable = find(is_usable)
  df = df[usable, symbols_to_keep]
  set_star_table(df, usable)
#end
  return df
end

function test_stellar_table() # TODO: Write test for stellar catalog functions
end
