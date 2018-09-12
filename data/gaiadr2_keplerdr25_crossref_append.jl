using FITSIO
using DataFrames
using CSV

crossref_tab = FITS("kepler_dr2_1arcsec.fits")[2];

#gaia_full = DataFrames.readtable("/sagan1/kepler/gaia_kepler_dr25_noblanklines.csv",separator=';',skipstart=243,nastrings=["NOT_AVAILABLE"]);
#deleterows!(gaia_full,1:2)

gaia_full = CSV.read("/sagan1/kepler/gaia_kepler_dr25_noblanklines.csv",delim=';',header=244,datarow=247,missingstring="NOT_AVAILABLE", rows_for_type_detect=1000)

(designation, source_id, gaia_ref_epoch, ra, ra_error, dec, dec_error, parallax, parallax_error, parallax_over_error, pmra, pmra_error, pmdec, pmdec_error, ra_dec_corr, ra_parallax_corr, ra_pmra_corr, ra_pmdec_corr, dec_parallax_corr, dec_pmra_corr, dec_pmdec_corr, parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, astrometric_chi2_al, astrometric_excess_noise, astrometric_excess_noise_sig, astrometric_primary_flag, duplicated_source, phot_g_mean_flux, phot_g_mean_flux_error, phot_g_mean_mag, phot_bp_mean_flux, phot_bp_mean_flux_error, phot_bp_mean_mag, phot_rp_mean_flux, phot_rp_mean_flux_error, phot_rp_mean_mag, bp_rp, bp_g, g_rp, radial_velocity, radial_velocity_error, phot_variable_flag, l, b, ecl_lon, ecl_lat, teff_val, teff_percentile_lower, teff_percentile_upper, a_g_val, a_g_percentile_lower, a_g_percentile_upper, e_bp_min_rp_val, e_bp_min_rp_percentile_lower, e_bp_min_rp_percentile_upper, radius_val, radius_percentile_lower, radius_percentile_upper, lum_val, lum_percentile_lower, lum_percentile_upper, r_est, r_lo, r_hi, r_length_prior, r_result_flag, r_modality_flag, kepid, tm_designation, ra_kepler, dec_kepler, kepmag, teff, teff_err1, teff_err2, teff_prov, logg, logg_err1, logg_err2, logg_prov, feh, feh_err1, feh_err2, feh_prov, radius, radius_err1, radius_err2, mass, mass_err1, mass_err2, prov_sec, nconfp, nkoi, ntce, jmag, hmag, kmag, planet_exist, kepler_gaia_ang_dist) = map(x->read(crossref_tab,x), ["designation", "source_id", "gaia_ref_epoch", "ra", "ra_error", "dec", "dec_error", "parallax", "parallax_error", "parallax_over_error", "pmra", "pmra_error", "pmdec", "pmdec_error", "ra_dec_corr", "ra_parallax_corr", "ra_pmra_corr", "ra_pmdec_corr", "dec_parallax_corr", "dec_pmra_corr", "dec_pmdec_corr", "parallax_pmra_corr", "parallax_pmdec_corr", "pmra_pmdec_corr", "astrometric_chi2_al", "astrometric_excess_noise", "astrometric_excess_noise_sig", "astrometric_primary_flag", "duplicated_source", "phot_g_mean_flux", "phot_g_mean_flux_error", "phot_g_mean_mag", "phot_bp_mean_flux", "phot_bp_mean_flux_error", "phot_bp_mean_mag", "phot_rp_mean_flux", "phot_rp_mean_flux_error", "phot_rp_mean_mag", "bp_rp", "bp_g", "g_rp", "radial_velocity", "radial_velocity_error", "phot_variable_flag", "l", "b", "ecl_lon", "ecl_lat", "teff_val", "teff_percentile_lower", "teff_percentile_upper", "a_g_val", "a_g_percentile_lower", "a_g_percentile_upper", "e_bp_min_rp_val", "e_bp_min_rp_percentile_lower", "e_bp_min_rp_percentile_upper", "radius_val", "radius_percentile_lower", "radius_percentile_upper", "lum_val", "lum_percentile_lower", "lum_percentile_upper", "r_est", "r_lo", "r_hi", "r_length_prior", "r_result_flag", "r_modality_flag", "kepid", "tm_designation", "ra_kepler", "dec_kepler", "kepmag", "teff", "teff_err1", "teff_err2", "teff_prov", "logg", "logg_err1", "logg_err2", "logg_prov", "feh", "feh_err1", "feh_err2", "feh_prov", "radius", "radius_err1", "radius_err2", "mass", "mass_err1", "mass_err2", "prov_sec", "nconfp", "nkoi", "ntce", "jmag", "hmag", "kmag", "planet?", "kepler_gaia_ang_dist"])

crossref_df = DataFrame(designation=designation, source_id=source_id, gaia_ref_epoch=gaia_ref_epoch, ra=ra, ra_error=ra_error, dec=dec, dec_error=dec_error, parallax=parallax, parallax_error=parallax_error, parallax_over_error=parallax_over_error, pmra=pmra, pmra_error=pmra_error, pmdec=pmdec, pmdec_error=pmdec_error, ra_dec_corr=ra_dec_corr, ra_parallax_corr=ra_parallax_corr, ra_pmra_corr=ra_pmra_corr, ra_pmdec_corr=ra_pmdec_corr, dec_parallax_corr=dec_parallax_corr, dec_pmra_corr=dec_pmra_corr, dec_pmdec_corr=dec_pmdec_corr, parallax_pmra_corr=parallax_pmra_corr, parallax_pmdec_corr=parallax_pmdec_corr, pmra_pmdec_corr=pmra_pmdec_corr)

dr2name = collect(skipmissing(gaia_full[:Source]))
astrometric_gof_al = collect(skipmissing(gaia_full[:gofAL]))
priam_flags = gaia_full[:fPriam]

crossref_ind = map(x->findfirst(dr2name,x),crossref_df[:source_id])

tmp_gof = Vector{Float64}(length(crossref_ind))
tmp_priam = Vector{Int64}(length(crossref_ind))

for i in 1:length(crossref_ind)
    tmp_gof[i] = astrometric_gof_al[crossref_ind[i]]
    if typeof(priam_flags[crossref_ind[i]]) == Missings.Missing
        tmp_priam[i] = -999
    else
        tmp_priam[i] = priam_flags[crossref_ind[i]]
    end
end

crossref_df[:astrometric_gof_al] = tmp_gof
#crossref_df[:astrometric_gof_al] = recode(crossref_df[:astrometric_gof_al], -999.0=>missing)

new_col = ["astrometric_chi2_al", "astrometric_excess_noise", "astrometric_excess_noise_sig", "astrometric_primary_flag", "duplicated_source", "phot_g_mean_flux", "phot_g_mean_flux_error", "phot_g_mean_mag", "phot_bp_mean_flux", "phot_bp_mean_flux_error", "phot_bp_mean_mag", "phot_rp_mean_flux", "phot_rp_mean_flux_error", "phot_rp_mean_mag", "bp_rp", "bp_g", "g_rp", "radial_velocity", "radial_velocity_error", "phot_variable_flag", "l", "b", "ecl_lon", "ecl_lat"]

tmp_col_names = vcat(names(crossref_df),[Symbol(i) for i in new_col])

crossref_df = hcat(crossref_df, astrometric_chi2_al, astrometric_excess_noise, astrometric_excess_noise_sig, astrometric_primary_flag, duplicated_source, phot_g_mean_flux, phot_g_mean_flux_error, phot_g_mean_mag, phot_bp_mean_flux, phot_bp_mean_flux_error, phot_bp_mean_mag, phot_rp_mean_flux, phot_rp_mean_flux_error, phot_rp_mean_mag, bp_rp, bp_g, g_rp, radial_velocity, radial_velocity_error, phot_variable_flag, l, b, ecl_lon, ecl_lat)

names!(crossref_df, tmp_col_names)

crossref_df[:priam_flags] = tmp_priam
crossref_df[:priam_flags] = recode(crossref_df[:priam_flags], -999=>missing)

new_col = ["teff_val", "teff_percentile_lower", "teff_percentile_upper", "a_g_val", "a_g_percentile_lower", "a_g_percentile_upper", "e_bp_min_rp_val", "e_bp_min_rp_percentile_lower", "e_bp_min_rp_percentile_upper", "radius_val", "radius_percentile_lower", "radius_percentile_upper", "lum_val", "lum_percentile_lower", "lum_percentile_upper", "r_est", "r_lo", "r_hi", "r_length_prior", "r_result_flag", "r_modality_flag", "kepid", "tm_designation", "ra_kepler", "dec_kepler", "kepmag", "teff", "teff_err1", "teff_err2", "teff_prov", "logg", "logg_err1", "logg_err2", "logg_prov", "feh", "feh_err1", "feh_err2", "feh_prov", "radius", "radius_err1", "radius_err2", "mass", "mass_err1", "mass_err2", "prov_sec", "nconfp", "nkoi", "ntce", "jmag", "hmag", "kmag", "planet_cand", "kepler_gaia_ang_dist"]

tmp_col_names = vcat(names(crossref_df),[Symbol(i) for i in new_col])

crossref_df = hcat(crossref_df, teff_val, teff_percentile_lower, teff_percentile_upper, a_g_val, a_g_percentile_lower, a_g_percentile_upper, e_bp_min_rp_val, e_bp_min_rp_percentile_lower, e_bp_min_rp_percentile_upper, radius_val, radius_percentile_lower, radius_percentile_upper, lum_val, lum_percentile_lower, lum_percentile_upper, r_est, r_lo, r_hi, r_length_prior, r_result_flag, r_modality_flag, kepid, tm_designation, ra_kepler, dec_kepler, kepmag, teff, teff_err1, teff_err2, teff_prov, logg, logg_err1, logg_err2, logg_prov, feh, feh_err1, feh_err2, feh_prov, radius, radius_err1, radius_err2, mass, mass_err1, mass_err2, prov_sec, nconfp, nkoi, ntce, jmag, hmag, kmag, planet_exist, kepler_gaia_ang_dist)

names!(crossref_df, tmp_col_names)

CSV.write("gaiadr2_keplerdr25_crossref.csv", crossref_df)
