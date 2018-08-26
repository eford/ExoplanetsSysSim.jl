using HDF5, JLD, DataFrames, CSV
using StatsBase, Polynomials, CurveFit

kep_filename = joinpath(Pkg.dir("ExoplanetsSysSim"), "data", "q1_q17_dr25_stellar.csv")
gaia_filename = joinpath(Pkg.dir("ExoplanetsSysSim"), "data", "gaiadr2_keplerdr25_crossref.csv")
stellar_catalog_file_out = joinpath(Pkg.dir("ExoplanetsSysSim"), "data", "q1q17_dr25_gaia_fgk.jld")

kep_df = CSV.read(kep_filename)
gaia_df = CSV.read(gaia_filename, rows_for_type_detect=30000)

dup_gaiaid = find(nonunique(DataFrame(x = gaia_df[:source_id])))
deleterows!(gaia_df, dup_gaiaid)

mag_diff = gaia_df[:phot_g_mean_mag].-gaia_df[:kepmag]
quant_arr = quantile(mag_diff, [0.067,0.933])   # 1.5-sigma cut
mag_match = find(x->quant_arr[1]<=x<=quant_arr[2], mag_diff)
gaia_df = gaia_df[mag_match,:]

gaia_col = [:kepid,:source_id,:parallax,:parallax_error,:astrometric_gof_al,:astrometric_excess_noise_sig,:bp_rp,:priam_flags,:teff_val,:teff_percentile_lower,:teff_percentile_upper,:radius_val,:radius_percentile_lower,:radius_percentile_upper,:lum_val,:lum_percentile_lower,:lum_percentile_upper]
df = join(kep_df, gaia_df[:,gaia_col], on=:kepid)
kep_df = nothing
gaia_df = nothing

df[:teff] = df[:teff_val]
df[:teff_err1] = df[:teff_percentile_upper].-df[:teff_val]
df[:teff_err2] = df[:teff_percentile_lower].-df[:teff_val]
delete!(df, :teff_val)
delete!(df, :teff_percentile_upper)
delete!(df, :teff_percentile_lower)
df[:radius_err1] = df[:radius_percentile_upper].-df[:radius_val]
df[:radius_err2] = df[:radius_percentile_lower].-df[:radius_val]
df[:radius] = df[:radius_val]
delete!(df, :radius_val)
delete!(df, :radius_percentile_upper)
delete!(df, :radius_percentile_lower)

not_binary = (df[:astrometric_gof_al] .<= 20) .& (df[:astrometric_excess_noise_sig] .<= 5)
astro_good = []
for x in 1:length(df[:kepid])
    if !(ismissing(df[x,:priam_flags]))
        pflag = string(df[x,:priam_flags])
         if (pflag[2] == '0') & (pflag[3] == '0')
             push!(astro_good, true)
         else
             push!(astro_good, false)
         end
     else
         push!(astro_good, false)
     end
end
astro_good = astro_good .& (df[:parallax_error] .< 0.05*df[:parallax])
planet_search = df[:kepmag] .<= 16.

has_mass = .! (ismissing.(df[:mass]) .| ismissing.(df[:mass_err1]) .| ismissing.(df[:mass_err2]))
has_radius = .! (ismissing.(df[:radius]) .| ismissing.(df[:radius_err1]) .| ismissing.(df[:radius_err2]))
has_dens = .! (ismissing.(df[:dens]) .| ismissing.(df[:dens_err1]) .| ismissing.(df[:dens_err2]))
has_rest = .! (ismissing.(df[:rrmscdpp04p5]) .| ismissing.(df[:dataspan]) .| ismissing.(df[:dutycycle]))

is_usable = has_radius .& has_mass .& has_rest .& has_dens .& astro_good .& not_binary .& planet_search

df = df[find(is_usable),:]

fgk_color = (0.5 .<= df[:bp_rp] .<= 1.7)
df = df[find(fgk_color),:]
is_FGK = []
coeff = [2.5,-3.6,0.9]
for i in 1:6
    is_FGK = (log10.(df[:lum_val]).< map(x->polyval(Poly(coeff),x),df[:bp_rp]) + log10(1.75))
    coeff = poly_fit(df[is_FGK,:bp_rp],log10.(df[is_FGK,:lum_val]),2)
end

# See options at: http://exoplanetarchive.ipac.caltech.edu/docs/API_keplerstellar_columns.html
# TODO SCI DETAIL or IMPORTANT?: Read in all CDPP's, so can interpolate?
symbols_to_keep = [ :kepid, :source_id, :mass, :mass_err1, :mass_err2, :radius, :radius_err1, :radius_err2, :dens, :dens_err1, :dens_err2, :rrmscdpp04p5, :dataspan, :dutycycle ]
delete!(df, [~(x in symbols_to_keep) for x in names(df)])    # delete columns that we won't be using anyway
FGK = find(is_FGK)
df = df[FGK, symbols_to_keep]
tmp_df = DataFrame()    
for col in names(df)
    tmp_df[col] = collect(skipmissing(df[col]))
end
df = tmp_df

save(stellar_catalog_file_out,"stellar_catalog", df)
