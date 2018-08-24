dir_list = [
"0_5-1_25_050-075_bin",
"0_5-1_25_075-100_bin",
"0_5-1_25_100-125_bin",
"0_5-1_25_125-150_bin",
"0_5-1_25_150-175_bin",
"0_5-1_25_175-200_bin",
"0_5-1_25_2-2_5_bin",
"0_5-1_25_2_5-3_bin",
"0_5-1_25_3-4_bin",
"0_5-1_25_4-6_bin",
"0_5-1_25_6-8_bin",
"0_5-1_25_8-12_bin",
"0_5-1_25_12-16_bin",
"1_25-2_5_050-075_bin",
"1_25-2_5_075-100_bin",
"1_25-2_5_100-125_bin",
"1_25-2_5_125-150_bin",
"1_25-2_5_150-175_bin",
"1_25-2_5_175-200_bin",
"1_25-2_5_2-2_5_bin",
"1_25-2_5_2_5-3_bin",
"1_25-2_5_3-4_bin",
"1_25-2_5_4-6_bin",
"1_25-2_5_6-8_bin",
"1_25-2_5_8-12_bin",
"1_25-2_5_12-16_bin",
"2_5-5_050-075_bin",
"2_5-5_075-100_bin",
"2_5-5_100-125_bin",
"2_5-5_125-150_bin",
"2_5-5_150-175_bin",
"2_5-5_175-200_bin",
"2_5-5_2-2_5_bin",
"2_5-5_2_5-3_bin",
"2_5-5_3-4_bin",
"2_5-5_4-6_bin",
"2_5-5_6-8_bin",
"2_5-5_8-12_bin",
"2_5-5_12-16_bin",
"5-10_050-075_bin",
"5-10_075-100_bin",
"5-10_100-125_bin",
"5-10_125-150_bin",
"5-10_150-175_bin",
"5-10_175-200_bin",
"5-10_2-2_5_bin",
"5-10_2_5-3_bin",
"5-10_3-4_bin",
"5-10_4-6_bin",
"5-10_6-8_bin",
"5-10_8-12_bin",
"5-10_12-16_bin",
"10-20_050-075_bin",
"10-20_075-100_bin",
"10-20_100-125_bin",
"10-20_125-150_bin",
"10-20_150-175_bin",
"10-20_175-200_bin",
"10-20_2-2_5_bin",
"10-20_2_5-3_bin",
"10-20_3-4_bin",
"10-20_4-6_bin",
"10-20_6-8_bin",
"10-20_8-12_bin",
"10-20_12-16_bin",
"20-40_050-075_bin",
"20-40_075-100_bin",
"20-40_100-125_bin",
"20-40_125-150_bin",
"20-40_150-175_bin",
"20-40_175-200_bin",
"20-40_2-2_5_bin",
"20-40_2_5-3_bin",
"20-40_3-4_bin",
"20-40_4-6_bin",
"20-40_6-8_bin",
"20-40_8-12_bin",
"20-40_12-16_bin",
"40-80_050-075_bin",
"40-80_075-100_bin",
"40-80_100-125_bin",
"40-80_125-150_bin",
"40-80_150-175_bin",
"40-80_175-200_bin",
"40-80_2-2_5_bin",
"40-80_2_5-3_bin",
"40-80_3-4_bin",
"40-80_4-6_bin",
"40-80_6-8_bin",
"40-80_8-12_bin",
"40-80_12-16_bin",
"80-160_050-075_bin",
"80-160_075-100_bin",
"80-160_100-125_bin",
"80-160_125-150_bin",
"80-160_150-175_bin",
"80-160_175-200_bin",
"80-160_2-2_5_bin",
"80-160_2_5-3_bin",
"80-160_3-4_bin",
"80-160_4-6_bin",
"80-160_6-8_bin",
"80-160_8-12_bin",
"80-160_12-16_bin",
"160-320_050-075_bin",
"160-320_075-100_bin",
"160-320_100-125_bin",
"160-320_125-150_bin",
"160-320_150-175_bin",
"160-320_175-200_bin",
"160-320_2-2_5_bin",
"160-320_2_5-3_bin",
"160-320_3-4_bin",
"160-320_4-6_bin",
"160-320_6-8_bin",
"160-320_8-12_bin",
"160-320_12-16_bin"
]

p_lim = [0.5, 1.25, 2.5, 5., 10., 20., 40., 80., 160., 320.]
r_lim = [0.5, 0.75, 1., 1.25, 1.5, 1.75, 2., 2.5, 3., 4., 6., 8., 12., 16.]

file_ind = 1
for p_ind = 1:(length(p_lim)-1)
    for r_ind = 1:(length(r_lim)-1)
        if !isdir(dir_list[file_ind])
            mkdir(dir_list[file_ind])
        end
        cd(dir_list[file_ind])
        f_param = open("param.in","w")
        write(f_param,"""stellar_catalog = "q1q17_dr25_gaia_fgk.jld"\n""")
        write(f_param,"""koi_catalog = "q1_q17_dr25_koi.csv"\n""")
        write(f_param,"num_targ_sim = 5000\n")
        write(f_param,string("p_bin_lim = [",p_lim[p_ind],", ",p_lim[p_ind+1],"]\n"))
        write(f_param,string("r_bin_lim = [",r_lim[r_ind],", ",r_lim[r_ind+1],"]\n"))
        write(f_param,"#rate_init = 1.0")
        file_ind += 1
        cd("..")
    end
end