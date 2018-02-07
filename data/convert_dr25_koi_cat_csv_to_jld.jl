using ExoplanetsSysSim
using HDF5, JLD, DataFrames

koi_catalog_file_in = joinpath(Pkg.dir("ExoplanetsSysSim"), "data", "q1_q17_dr25_koi.csv")
koi_catalog_file_out = replace(koi_catalog_file_in,r".csv$",".jld")

(csv_data,csv_header) =  readcsv(koi_catalog_file_in,header=true)
# Lookup header columns, since DataFrames doesn't like this file
kepid_idx = findfirst(x->x=="kepid",csv_header)
koi_period_idx = findfirst(x->x=="koi_period",csv_header)
koi_time0bk_idx = findfirst(x->x=="koi_time0bk",csv_header)
koi_depth_idx = findfirst(x->x=="koi_depth",csv_header)
koi_duration_idx = findfirst(x->x=="koi_duration",csv_header)
koi_ror_idx = findfirst(x->x=="koi_ror",csv_header)
koi_period_err1_idx = findfirst(x->x=="koi_period_err1",csv_header)
koi_time0bk_err1_idx = findfirst(x->x=="koi_time0bk_err1",csv_header)
koi_depth_err1_idx = findfirst(x->x=="koi_depth_err1",csv_header)
koi_duration_err1_idx = findfirst(x->x=="koi_duration_err1",csv_header)
koi_period_err2_idx = findfirst(x->x=="koi_period_err2",csv_header)
koi_time0bk_err2_idx = findfirst(x->x=="koi_time0bk_err2",csv_header)
koi_depth_err2_idx = findfirst(x->x=="koi_depth_err2",csv_header)
koi_duration_err2_idx = findfirst(x->x=="koi_duration_err2",csv_header)
koi_disposition_idx = findfirst(x->x=="koi_disposition",csv_header)
koi_pdisposition_idx = findfirst(x->x=="koi_pdisposition",csv_header)
symbols_to_keep = [kepid_idx,koi_period_idx,koi_time0bk_idx,koi_depth_idx,koi_duration_idx,koi_ror_idx,koi_period_err1_idx,koi_time0bk_err1_idx,koi_depth_err1_idx,koi_duration_err1_idx,koi_period_err2_idx,koi_time0bk_err2_idx,koi_depth_err2_idx,koi_duration_err2_idx]
# Choose which KOIs to keep
#is_cand = (csv_data[:,koi_disposition_idx] .== "CONFIRMED") | (csv_data[:,koi_disposition_idx] .== "CANDIDATE")
is_cand = (csv_data[:,koi_pdisposition_idx] .== "CANDIDATE")

idx_keep = is_cand & !isna(csv_data[:,koi_ror_idx]) & ([typeof(x) for x in csv_data[:,koi_ror_idx]] .== Float64)
idx_keep = idx_keep & !isna(csv_data[:,koi_period_err1_idx]) & ([typeof(x) for x in csv_data[:,koi_period_err1_idx]] .== Float64) # DR25 catalog missing uncertainties for some candidates
csv_data = csv_data[idx_keep,symbols_to_keep]

save(koi_catalog_file_out,"koi_catalog", csv_data)








 
