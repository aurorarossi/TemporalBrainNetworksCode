using NPZ, Random, Distributions, LinearAlgebra, Permutations, Plots, JLD2, Statistics, DelimitedFiles, Dierckx

# LOAD DATA 
function load_fMRI_data(subjects, thresholds)
    # initialize variables
    len_sub = length(subjects)
    len_thre = length(thresholds)
    average_deg_all = zeros(len_sub, len_thre)
    clustering_all = zeros(len_sub, len_thre)
    path_len_all = zeros(len_sub, len_thre)
    temp_clustering_all = zeros(len_sub, len_thre)
    # load data
    for i in 1:len_sub
        @load "/data/Hyperbrain/$(subjects[i])/$(subjects[i])_schaefer_300_thre_02_005_098_OR_RT_RE_ws60_wo30.jld2" av cl pl tcc
        average_deg_all[i, :] = av
        clustering_all[i, :] = cl
        path_len_all[i, :] = pl
        temp_clustering_all[i, :] = tcc
    end
    return average_deg_all, clustering_all, path_len_all, temp_clustering_all
end

function my_quantile(len_thre, variable)
    quantile_variable = zeros(3, len_thre)
    for j in 1:len_thre
        quantile_variable[:, j] = quantile(variable[:, j])[2:4]
    end
    return quantile_variable
end

# COMPUTE QUANTILES
function compute_quantiles(average_deg_all, clustering_all, path_len_all, temp_clustering_all, thresholds)
    len_thre = length(thresholds)
    # initialize variables
    average_deg_quantile = zeros(len_thre)
    for j in 1:len_thre
        average_deg_quantile[j] = mean(average_deg_all[:, j])
    end
    clustering_quantile = my_quantile(len_thre, clustering_all)
    path_len_quantile = my_quantile(len_thre, path_len_all)
    small_world_quantile = my_quantile(len_thre, clustering_all ./ path_len_all)
    small_world_SB_quantile = my_quantile(len_thre, temp_clustering_all ./ path_len_all)
    temp_clustering_quantile = my_quantile(len_thre, temp_clustering_all)
    return average_deg_quantile, clustering_quantile, path_len_quantile, small_world_quantile, small_world_SB_quantile, temp_clustering_quantile
end

function load_hyperbolic_data(len_R=length(Rrange))
    αrange = collect(0.5:0.025:1.2)
    velocities = collect(0.1:0.1:0.9)
    len_vel = length(velocities)
    len_α = length(αrange)
    average_deg_HY = zeros(len_vel, len_R, len_α)
    clustering_HY = zeros(len_vel, len_R, len_α)
    path_len_HY = zeros(len_vel, len_R, len_α)
    temp_clustering_HY = zeros(len_vel, len_R, len_α)
    small_world_HY = zeros(len_vel, len_R, len_α)
    small_world_SB_HY = zeros(len_vel, len_R, len_α)
    for a in 1:len_α
        @load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/hyperbolic_300_alpha_$(round(Int,αrange[a]*100))_R_05_05_14_v_01_01_09_correct.jld2" avHY clHY plHY tccHY
        average_deg_HY[:, :, a] = avHY[:, :]
        clustering_HY[:, :, a] = clHY[:, :]
        path_len_HY[:, :, a] = plHY[:, :]
        temp_clustering_HY[:, :, a] = tccHY[:, :]
        small_world_HY[:, :, a] = clHY[:, :] ./ plHY[:, :]
        small_world_SB_HY[:, :, a] = tccHY[:, :] ./ plHY[:, :]
    end
    return average_deg_HY, clustering_HY, path_len_HY, temp_clustering_HY, small_world_HY, small_world_SB_HY
end

function compute_integral(variable_fMRI, variable_HY, len_vel=length(velocities), len_α=length(αrange))
    integral = zeros(len_vel, len_α)


end


function main()
    thresholds = append!(collect(0.2:0.05:0.90), collect(0.92:0.02:0.98))
    subjects = readdlm("TemporalBrainNetworksCode/src/filtered-subjects-mod.txt", Int)
    Rrange=collect(0.5:0.5:14)
    average_deg_all, clustering_all, path_len_all, temp_clustering_all = load_fMRI_data(subjects, thresholds)
    average_deg_quantile, clustering_quantile, path_len_quantile, small_world_quantile, small_world_SB_quantile, temp_clustering_quantile = compute_quantiles(average_deg_all, clustering_all, path_len_all, temp_clustering_all, thresholds)
    average_deg_HY, clustering_HY, path_len_HY, temp_clustering_HY, small_world_HY, small_world_SB_HY = load_hyperbolic_data()

    p = plot(average_deg_quantile, small_world_quantile[2, :], label="Real data 1050 Schaefer", lw=3)
    plot!(average_deg_HY[1,:,1], small_world_HY[1, :,1], label="", lw=3, ls=:dash)
end