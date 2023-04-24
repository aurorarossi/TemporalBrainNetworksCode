# LOAD DATA 
function load_fMRI_data(subjects, thresholds)
    # initialize variables
    len_sub = length(subjects)
    len_thre = length(thresholds)
    average_deg_all = zeros(len_sub, len_thre)
    clustering_all = zeros(len_sub, len_thre)
    path_len_all = zeros(len_sub, len_thre)
    temp_corr_all = zeros(len_sub, len_thre)
    # load data
    for i in 1:len_sub
        @load "/data/Hyperbrain/$(subjects[i])/$(subjects[i])_schaefer_300_thre_02_005_098_OR_RT_RE_ws60_wo30.jld2" av cl pl tcc
        average_deg_all[i, :] = av
        clustering_all[i, :] = cl
        path_len_all[i, :] = pl
        temp_corr_all[i, :] = tcc
    end
    return average_deg_all, clustering_all, path_len_all, temp_corr_all
end

function load_fMRI_data_RT(subjects, thresholds)
    # initialize variables
    len_sub = length(subjects)
    len_thre = length(thresholds)
    average_deg_all = zeros(len_sub, len_thre)
    clustering_all = zeros(len_sub, len_thre)
    path_len_all = zeros(len_sub, len_thre)
    temp_corr_all = zeros(len_sub, len_thre)
    # load data
    for i in 1:len_sub
        @load "/data/Hyperbrain/$(subjects[i])/$(subjects[i])_schaefer_300_thre_02_005_098_OR_RT_RE_ws60_wo30.jld2" avRT clRT plRT tccRT
        average_deg_all[i, :] = avRT
        clustering_all[i, :] = clRT
        path_len_all[i, :] = plRT
        temp_corr_all[i, :] = tccRT
    end
    return average_deg_all, clustering_all, path_len_all, temp_corr_all
end

function load_fMRI_data_RE(subjects, thresholds)
    # initialize variables
    len_sub = length(subjects)
    len_thre = length(thresholds)
    average_deg_all = zeros(len_sub, len_thre)
    clustering_all = zeros(len_sub, len_thre)
    path_len_all = zeros(len_sub, len_thre)
    temp_corr_all = zeros(len_sub, len_thre)
    # load data
    for i in 1:len_sub
        @load "/data/Hyperbrain/$(subjects[i])/$(subjects[i])_schaefer_300_thre_02_005_098_OR_RT_RE_ws60_wo30.jld2" avRE clRE plRE tccRE
        average_deg_all[i, :] = avRE
        clustering_all[i, :] = clRE
        path_len_all[i, :] = plRE
        temp_corr_all[i, :] = tccRE
    end
    return average_deg_all, clustering_all, path_len_all, temp_corr_all
end

function load_euclidean_data(len)
    # initialize variables
    average_deg_all = zeros(10, len)
    clustering_all = zeros(10, len)
    path_len_all = zeros(10, len)
    temp_corr_all = zeros(10, len)
    # load data
    for i in 1:10
        @load "/user/aurossi/home/mri_networks/TemporalBrainNetworksCode/data_correct/euclidean_300_thre_02_005_098_v02_$(i).jld2" avEU clEU plEU tccEU
        average_deg_all[i, :] = avEU
        clustering_all[i, :] = clEU
        path_len_all[i, :] = plEU
        temp_corr_all[i, :] = tccEU
    end
    return average_deg_all, clustering_all, path_len_all, temp_corr_all
end
function load_euclideanBorder_data(len)
    # initialize variables
    average_deg_all = zeros(10, len)
    clustering_all = zeros(10, len)
    path_len_all = zeros(10, len)
    temp_corr_all = zeros(10, len)
    # load data
    for i in 1:10
        @load "/user/aurossi/home/mri_networks/TemporalBrainNetworksCode/data_correct/euclideanBorder_300_r_0_005_08_v01_$(i).jld2" avEU clEU plEU tccEU
        average_deg_all[i, :] = avEU
        clustering_all[i, :] = clEU
        path_len_all[i, :] = plEU
        temp_corr_all[i, :] = tccEU
    end
    return average_deg_all, clustering_all, path_len_all, temp_corr_all
end

function load_hyperbolic_data(len)
    # initialize variables
    average_deg_all = zeros(10, len)
    clustering_all = zeros(10, len)
    path_len_all = zeros(10, len)
    temp_corr_all = zeros(10, len)
    # load data
    for i in 1:10
        @load "/user/aurossi/home/mri_networks/TemporalBrainNetworksCode/data_correct/hyperbolic_300_alpha_65_R_0_05_18_v_60_$(i).jld2" avHY clHY plHY tccHY
        average_deg_all[i, :] = avHY
        clustering_all[i, :] = clHY
        path_len_all[i, :] = plHY
        temp_corr_all[i, :] = tccHY
    end
    return average_deg_all, clustering_all, path_len_all, temp_corr_all
end

function load_hyperbolic_SB_data(len)
    # initialize variables
    average_deg_all = zeros(10, len)
    clustering_all = zeros(10, len)
    path_len_all = zeros(10, len)
    temp_corr_all = zeros(10, len)
    # load data
    for i in 1:10
        @load "/user/aurossi/home/mri_networks/TemporalBrainNetworksCode/data_correct/hyperbolic_300_alpha_80_R_0_05_18_v_90_$(i).jld2" avHY clHY plHY tccHY
        average_deg_all[i, :] = avHY
        clustering_all[i, :] = clHY
        path_len_all[i, :] = plHY
        temp_corr_all[i, :] = tccHY
    end
    return average_deg_all, clustering_all, path_len_all, temp_corr_all
end

function my_quantile(len, variable)
    quantile_variable = zeros(3, len)
    for j in 1:len
        quantile_variable[:, j] = quantile(variable[:, j])[2:4]
    end
    return quantile_variable
end

function compute_quantiles(average_deg_all, clustering_all, path_len_all, temp_corr_all, len)
    # initialize variables
    average_deg_quantile = zeros(len)
    for j in 1:len
        average_deg_quantile[j] = mean(average_deg_all[:, j])
    end
    clustering_quantile = my_quantile(len, clustering_all)
    path_len_quantile = my_quantile(len, path_len_all)
    small_world_quantile = my_quantile(len, clustering_all ./ path_len_all)
    small_world_SB_quantile = my_quantile(len, temp_corr_all ./ path_len_all)
    temp_corr_quantile = my_quantile(len, temp_corr_all)
    return average_deg_quantile, clustering_quantile, path_len_quantile, small_world_quantile, small_world_SB_quantile, temp_corr_quantile
end

function load_and_quantiles()
    thresholds = append!(collect(0.2:0.05:0.90), collect(0.92:0.02:0.98))
    subjects = readdlm("src/filtered-subjects-mod.txt", Int)
    radiusRange = collect(0:0.05:0.8)
    Rrange = collect(0:0.5:18)
    average_deg_all, clustering_all, path_len_all, temp_corr_all = load_fMRI_data(subjects, thresholds)
    average_deg_all_RT, clustering_all_RT, path_len_all_RT, temp_corr_all_RT = load_fMRI_data_RT(subjects, thresholds)
    average_deg_all_RE, clustering_all_RE, path_len_all_RE, temp_corr_all_RE = load_fMRI_data_RE(subjects, thresholds)
    average_deg_all_EU, clustering_all_EU, path_len_all_EU, temp_corr_all_EU = load_euclidean_data(length(thresholds))
    average_deg_all_EUB, clustering_all_EUB, path_len_all_EUB, temp_corr_all_EUB = load_euclideanBorder_data(length(radiusRange))
    average_deg_all_HY, clustering_all_HY, path_len_all_HY, temp_corr_all_HY = load_hyperbolic_data(length(Rrange))
    average_deg_all_HYSB, clustering_all_HYSB, path_len_all_HYSB, temp_corr_all_HYSB = load_hyperbolic_SB_data(length(Rrange))
    average_deg_quantile, clustering_quantile, path_len_quantile, small_world_quantile, small_world_SB_quantile, temp_corr_quantile = compute_quantiles(average_deg_all, clustering_all, path_len_all, temp_corr_all, length(thresholds))
    average_deg_quantile_RT, clustering_quantile_RT, path_len_quantile_RT, small_world_quantile_RT, small_world_SB_quantile_RT, temp_corr_quantile_RT = compute_quantiles(average_deg_all_RT, clustering_all_RT, path_len_all_RT, temp_corr_all_RT, length(thresholds))
    average_deg_quantile_RE, clustering_quantile_RE, path_len_quantile_RE, small_world_quantile_RE, small_world_SB_quantile_RE, temp_corr_quantile_RE = compute_quantiles(average_deg_all_RE, clustering_all_RE, path_len_all_RE, temp_corr_all_RE, length(thresholds))
    average_deg_quantile_EU, clustering_quantile_EU, path_len_quantile_EU, small_world_quantile_EU, small_world_SB_quantile_EU, temp_corr_quantile_EU = compute_quantiles(average_deg_all_EU, clustering_all_EU, path_len_all_EU, temp_corr_all_EU, length(thresholds))
    average_deg_quantile_EUB, clustering_quantile_EUB, path_len_quantile_EUB, small_world_quantile_EUB, small_world_SB_quantile_EUB, temp_corr_quantile_EUB = compute_quantiles(average_deg_all_EUB, clustering_all_EUB, path_len_all_EUB, temp_corr_all_EUB, length(radiusRange))
    average_deg_quantile_HY, clustering_quantile_HY, path_len_quantile_HY, small_world_quantile_HY, small_world_SB_quantile_HY, temp_corr_quantile_HY = compute_quantiles(average_deg_all_HY, clustering_all_HY, path_len_all_HY, temp_corr_all_HY, length(Rrange))
    average_deg_quantile_HYSB, clustering_quantile_HYSB, path_len_quantile_HYSB, small_world_quantile_HYSB, small_world_SB_quantile_HYSB, temp_corr_quantile_HYSB = compute_quantiles(average_deg_all_HYSB, clustering_all_HYSB, path_len_all_HYSB, temp_corr_all_HYSB, length(Rrange))
    return average_deg_quantile, clustering_quantile, path_len_quantile, small_world_quantile, small_world_SB_quantile, temp_corr_quantile, average_deg_quantile_RT, clustering_quantile_RT, path_len_quantile_RT, small_world_quantile_RT, small_world_SB_quantile_RT, temp_corr_quantile_RT, average_deg_quantile_RE, clustering_quantile_RE, path_len_quantile_RE, small_world_quantile_RE, small_world_SB_quantile_RE, temp_corr_quantile_RE, average_deg_quantile_EU, clustering_quantile_EU, path_len_quantile_EU, small_world_quantile_EU, small_world_SB_quantile_EU, temp_corr_quantile_EU, average_deg_quantile_EUB, clustering_quantile_EUB, path_len_quantile_EUB, small_world_quantile_EUB, small_world_SB_quantile_EUB, temp_corr_quantile_EUB, average_deg_quantile_HY, clustering_quantile_HY, path_len_quantile_HY, small_world_quantile_HY, small_world_SB_quantile_HY, temp_corr_quantile_HY, average_deg_quantile_HYSB, clustering_quantile_HYSB, path_len_quantile_HYSB, small_world_quantile_HYSB, small_world_SB_quantile_HYSB, temp_corr_quantile_HYSB
end

function compute_integral(average_deg_MOD, measure_MOD, average_fMRI, measure_fMRI)
    range = (1.5:0.01:160)
    spline = linear_interpolation(sort(average_deg_MOD), measure_MOD[sortperm(average_deg_MOD)])
    spline_fMRI = linear_interpolation(sort(average_fMRI), measure_fMRI[sortperm(average_fMRI)])
    diffmy = abs.(spline(range) .- spline_fMRI(range))
    diffspline = linear_interpolation(range, diffmy)
    integral, _ = quadgk(diffspline, range[1], range[end])
    return integral
end

average_deg_quantile, clustering_quantile, path_len_quantile, small_world_quantile, small_world_SB_quantile, temp_corr_quantile, average_deg_quantile_RT, clustering_quantile_RT, path_len_quantile_RT, small_world_quantile_RT, small_world_SB_quantile_RT, temp_corr_quantile_RT, average_deg_quantile_RE, clustering_quantile_RE, path_len_quantile_RE, small_world_quantile_RE, small_world_SB_quantile_RE, temp_corr_quantile_RE, average_deg_quantile_EU, clustering_quantile_EU, path_len_quantile_EU, small_world_quantile_EU, small_world_SB_quantile_EU, temp_corr_quantile_EU, average_deg_quantile_EUB, clustering_quantile_EUB, path_len_quantile_EUB, small_world_quantile_EUB, small_world_SB_quantile_EUB, temp_corr_quantile_EUB, average_deg_quantile_HY, clustering_quantile_HY, path_len_quantile_HY, small_world_quantile_HY, small_world_SB_quantile_HY, temp_corr_quantile_HY, average_deg_quantile_HYSB, clustering_quantile_HYSB, path_len_quantile_HYSB, small_world_quantile_HYSB, small_world_SB_quantile_HYSB, temp_corr_quantile_HYSB = load_and_quantiles()

#small worldness SB
println("SWSB")
println(compute_integral(average_deg_quantile_HYSB, small_world_SB_quantile_HYSB[2,:], average_deg_quantile, small_world_SB_quantile[2,:]))
println(compute_integral(average_deg_quantile_HY, small_world_SB_quantile_HY[2,:], average_deg_quantile, small_world_SB_quantile[2,:]))
println(compute_integral(average_deg_quantile_EUB, small_world_SB_quantile_EUB[2,:], average_deg_quantile, small_world_SB_quantile[2,:]))
println(compute_integral(average_deg_quantile_EU, small_world_SB_quantile_EU[2,:], average_deg_quantile, small_world_SB_quantile[2,:]))
println(compute_integral(average_deg_quantile_RE, small_world_SB_quantile_RE[2,:], average_deg_quantile, small_world_SB_quantile[2,:]))
println(compute_integral(average_deg_quantile_RT, small_world_SB_quantile_RT[2,:], average_deg_quantile, small_world_SB_quantile[2,:]))

#small worldness
println("SW")
println("HYSB",compute_integral(average_deg_quantile_HYSB, small_world_quantile_HYSB[2,:], average_deg_quantile, small_world_quantile[2,:]))
println("HY",compute_integral(average_deg_quantile_HY, small_world_quantile_HY[2,:], average_deg_quantile, small_world_quantile[2,:]))
println("EUB",compute_integral(average_deg_quantile_EUB, small_world_quantile_EUB[2,:], average_deg_quantile, small_world_quantile[2,:]))
println("EU",compute_integral(average_deg_quantile_EU, small_world_quantile_EU[2,:], average_deg_quantile, small_world_quantile[2,:]))
println("RE",compute_integral(average_deg_quantile_RE, small_world_quantile_RE[2,:], average_deg_quantile, small_world_quantile[2,:]))
println("RT",compute_integral(average_deg_quantile_RT, small_world_quantile_RT[2,:], average_deg_quantile, small_world_quantile[2,:]))

# #clustering
# compute_integral(average_deg_quantile_HYSB, clustering_quantile_HYSB, average_deg_quantile, clustering_quantile)
# compute_integral(average_deg_quantile_HY, clustering_quantile_HY, average_deg_quantile, clustering_quantile)
# compute_integral(average_deg_quantile_EUB, clustering_quantile_EUB, average_deg_quantile, clustering_quantile)
# compute_integral(average_deg_quantile_EU, clustering_quantile_EU, average_deg_quantile, clustering_quantile)
# compute_integral(average_deg_quantile_RE, clustering_quantile_RE, average_deg_quantile, clustering_quantile)
# compute_integral(average_deg_quantile_RT, clustering_quantile_RT, average_deg_quantile, clustering_quantile)

# #temporal correlation
# compute_integral(average_deg_quantile_HYSB, temp_corr_quantile_HYSB, average_deg_quantile, temp_corr_quantile)
# compute_integral(average_deg_quantile_HY, temp_corr_quantile_HY, average_deg_quantile, temp_corr_quantile)
# compute_integral(average_deg_quantile_EUB, temp_corr_quantile_EUB, average_deg_quantile, temp_corr_quantile)
# compute_integral(average_deg_quantile_EU, temp_corr_quantile_EU, average_deg_quantile, temp_corr_quantile)
# compute_integral(average_deg_quantile_RE, temp_corr_quantile_RE, average_deg_quantile, temp_corr_quantile)
# compute_integral(average_deg_quantile_RT, temp_corr_quantile_RT, average_deg_quantile, temp_corr_quantile)

# #path length
# compute_integral(average_deg_quantile_HYSB, path_len_quantile_HYSB, average_deg_quantile, path_len_quantile)
# compute_integral(average_deg_quantile_HY, path_len_quantile_HY, average_deg_quantile, path_len_quantile)
# compute_integral(average_deg_quantile_EUB, path_len_quantile_EUB, average_deg_quantile, path_len_quantile)
# compute_integral(average_deg_quantile_EU, path_len_quantile_EU, average_deg_quantile, path_len_quantile)
# compute_integral(average_deg_quantile_RE, path_len_quantile_RE, average_deg_quantile, path_len_quantile)
# compute_integral(average_deg_quantile_RT, path_len_quantile_RT, average_deg_quantile, path_len_quantile)

