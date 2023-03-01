using NPZ, Random, Distributions, LinearAlgebra, Permutations, Plots, JLD2, Statistics, DelimitedFiles, Dierckx

thresholds = append!(collect(0.2:0.05:0.90), collect(0.92:0.02:0.98))
subjects = readdlm("TemporalBrainNetworksCode/src/filtered-subjects-mod.txt", Int)

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

function main()
    average_deg_all, clustering_all, path_len_all, temp_clustering_all = load_fMRI_data(subjects, thresholds)
    average_deg_quantile, clustering_quantile, path_len_quantile, small_world_quantile, small_world_SB_quantile, temp_clustering_quantile = compute_quantiles(average_deg_all, clustering_all, path_len_all, temp_clustering_all, thresholds)
    p=plot(average_deg_quantile,small_world_quantile[2,:],label="Real data 1050 Schaefer",lw=3)
end
