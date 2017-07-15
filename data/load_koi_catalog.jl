using DataFrames
using ExoplanetsSysSim

function make_actual_stellar_catalog()
  stellar_catalog_file_in = ascii(joinpath(Pkg.dir("ExoplanetsSysSim"), "data", "q1_q17_dr25_stellar.csv"))
  df =  readtable(stellar_catalog_file_in)
  has_mass = ! (isna(df[:mass]) | isna(df[:mass_err1]) | isna(df[:mass_err2]))
  has_radius = ! (isna(df[:radius]) | isna(df[:radius_err1]) | isna(df[:radius_err2]))
  has_teff = ! (isna(df[:teff]) | isna(df[:teff_err1]) | isna(df[:teff_err2]))
  has_logg = ! (isna(df[:logg]) | isna(df[:logg_err1]) | isna(df[:logg_err2]))
  # TODO:  Add additional filters to match stellar catalog of Christian et al
  is_usable = has_mass & has_radius & has_teff & has_logg
  symbols_to_keep = [ :kepid, :mass, :radius, :teff, :logg ]
  df = df[find(is_usable), symbols_to_keep]
  return df
end

# df_star is assumed to have fields kepid, mass and radius for all targets in the survey
function make_actual_planet_candidate_catalog(df_star::DataFrame)
  koi_catalog_file_in = ascii(joinpath(Pkg.dir("ExoplanetsSysSim"), "data", "q1_q17_dr25_koi.csv"))
  (csv_data,csv_header) =  readcsv(koi_catalog_file_in,header=true)
  # Lookup header columns, since DataFrames doesn't like this file
  kepid_idx = findfirst(x->x=="kepid",csv_header)
  koi_period_idx = findfirst(x->x=="koi_period",csv_header)
  koi_time0bk_idx = findfirst(x->x=="koi_time0bk",csv_header)
  koi_depth_idx = findfirst(x->x=="koi_depth",csv_header)
  koi_duration_idx = findfirst(x->x=="koi_duration",csv_header)
  koi_period_err1_idx = findfirst(x->x=="koi_period_err1",csv_header)
  koi_time0bk_err1_idx = findfirst(x->x=="koi_time0bk_err1",csv_header)
  koi_depth_err1_idx = findfirst(x->x=="koi_depth_err1",csv_header)
  koi_duration_err1_idx = findfirst(x->x=="koi_duration_err1",csv_header)
  koi_period_err2_idx = findfirst(x->x=="koi_period_err2",csv_header)
  koi_time0bk_err2_idx = findfirst(x->x=="koi_time0bk_err2",csv_header)
  koi_depth_err2_idx = findfirst(x->x=="koi_depth_err2",csv_header)
  koi_duration_err2_idx = findfirst(x->x=="koi_duration_err2",csv_header)
  koi_disposition_idx = findfirst(x->x=="koi_disposition",csv_header)
  # Choose which KOIs to keep
  idx_keep = (csv_data[:,koi_disposition_idx] .== "CANDIDATE") | (csv_data[:,koi_disposition_idx] .== "CONFIRMED" ) # TODO:  Before publication, check that this is the best disposition field to use
  # TODO:  Add additional filters to match what planet candidates are to be kept Christian et al
  csv_data = csv_data[idx_keep,:]
 
  output = KeplerObsCatalog([])

  for (j,kepid) in enumerate(df_star[:kepid])
     plids = find(x->x==kepid,csv_data[:,kepid_idx])
	 num_pl = length(plids)
	 if num_pl ==0 continue end
	 target_obs = KeplerTargetObs(num_pl)
	 target_obs.star = ExoplanetsSysSim.StarObs(df_star[j,:radius],df_star[j,:mass])
	 for (plid,i) in enumerate(plids)
	   target_obs.obs[plid] = ExoplanetsSysSim.TransitPlanetObs(csv_data[i,koi_period_idx],csv_data[i,koi_time0bk_idx],csv_data[i,koi_depth_idx],csv_data[i,koi_duration_idx])
	   target_obs.sigma[plid] = ExoplanetsSysSim.TransitPlanetObs((abs(csv_data[i,koi_period_err1_idx])+abs(csv_data[i,koi_period_err2_idx]))/2,(abs(csv_data[i,koi_time0bk_err1_idx])+abs(csv_data[i,koi_time0bk_err2_idx]))/2,(abs(csv_data[i,koi_depth_err1_idx])+abs(csv_data[i,koi_depth_err2_idx]))/2,(abs(csv_data[i,koi_duration_err1_idx])+abs(csv_data[i,koi_duration_err2_idx]))/2)
	   target_obs.prob_detect = ExoplanetsSysSim.OneObserverSystemDetectionProbs( ones(num_pl), ones(num_pl,num_pl), ones(num_pl), fill(Array(Int64,0), 1) )
# ExoplanetsSysSim.OneObserverSystemDetectionProbs(num_pl)
	   push!(output.target,target_obs)
	 end	 
  end
  return output
end

# Demo of calling the two functions
#df_star = make_actual_stellar_catalog()
#planet_candidate_catalog = make_actual_planet_candidate_catalog(df_star)

