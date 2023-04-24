using NPZ, Random, Distributions, LinearAlgebra, Permutations, Plots, JLD2, Statistics, DelimitedFiles, Interpolations, QuadGK

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

# INITIALIZE AND COMPUTE QUANTILES
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

# LOAD DATA HYPERBOLIC
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

#COMPUTE INTEGRAL DISTANCE HYPERBOLIC AND DATA CURVES
function compute_integral_hyperbolic(average_deg_quantile, average_deg_HY, variable_fMRI, variable_HY, len_vel=length(velocities), len_α=length(αrange))
    integral = zeros(len_vel, len_α)
    range = collect(10:0.01:160)
    for a in 1:len_α
        for i in 1:len_vel
            spline = linear_interpolation(sort(average_deg_HY[i, :, a]), variable_HY[i, sortperm(average_deg_HY[i, :, a]), a])
            spline_fMRI = linear_interpolation(sort(average_deg_quantile), variable_fMRI[2, sortperm(average_deg_quantile)])
            diffmy = abs.(spline(range) .- spline_fMRI(range))
            diffspline = linear_interpolation(range, diffmy)
            integral[i, a], _ = quadgk(diffspline, range[1], range[end])
        end
    end
    return integral
end

function load_euclidean_data()
    @load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/euclidean_300_thre_02_005_098.jld2" avEU clEU plEU tccEU
    return avEU, clEU, plEU, tccEU
end

function load_euclideanBorder_data()
    @load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/euclideanBorder_300_r_0_005_08.jld2" avEU clEU plEU tccEU
    return avEU, clEU, plEU, tccEU
end

function compute_integral_euclidean(average_deg_quantile, variable_fMRI, average_deg_EU, variable_EU, len_vel, len_thre)
    integral_euclidean = zeros(len_vel)
    range = collect(10:0.01:160)
    for i in 1:len_vel
        spline = linear_interpolation(sort(average_deg_EU[i, :]), variable_EU[i, sortperm(average_deg_EU[i, :])])
        spline_fMRI = linear_interpolation(sort(average_deg_quantile), variable_fMRI[2, sortperm(average_deg_quantile)])
        diffmy = abs.(spline(range) .- spline_fMRI(range))
        diffspline = linear_interpolation(range, diffmy)
        integral_euclidean[i], _ = quadgk(diffspline, range[1], range[end])
    end
    return integral_euclidean
end

function find_minimum_v(integral, velocities)
    minimum_area = minimum(integral)
    indices = argmin(integral)
    best_velocity = velocities[indices]
    return minimum_area, best_velocity, indices
end

function find_minima_α_v(integral, αrange=collect(0.5:0.025:1.2), velocities=collect(0.1:0.1:0.9))
    minimum_area = minimum(integral)
    indices = argmin(integral)
    best_velocity = velocities[indices[1]]
    best_α = αrange[indices[2]]
    return minimum_area, best_velocity, best_α, indices
end



function main()
    thresholds = append!(collect(0.2:0.05:0.90), collect(0.92:0.02:0.98))
    subjects = readdlm("src/filtered-subjects-mod.txt", Int)
    Rrange = collect(0.5:0.5:14)
    average_deg_all, clustering_all, path_len_all, temp_clustering_all = load_fMRI_data(subjects, thresholds)
    average_deg_quantile, clustering_quantile, path_len_quantile, small_world_quantile, small_world_SB_quantile, temp_clustering_quantile = compute_quantiles(average_deg_all, clustering_all, path_len_all, temp_clustering_all, thresholds)
    len_R = length(Rrange)
    velocities = collect(0.1:0.1:0.9)
    αrange = collect(0.5:0.025:1.2)
    average_deg_HY, clustering_HY, path_len_HY, temp_clustering_HY, small_world_HY, small_world_SB_HY = load_hyperbolic_data(len_R)
    integral = compute_integral_hyperbolic(average_deg_quantile, average_deg_HY, small_world_SB_quantile, small_world_SB_HY, length(velocities), length(αrange))
    minimum_area, best_velocity, best_α, indices = find_minima_α_v(integral, αrange, velocities)
    println("Minimum area: ", minimum_area)
    println("Best velocity: ", best_velocity)
    println("Best α: ", best_α)
    println("Indices: ", indices)
    avEU, clEU, plEU, tccEU = load_euclidean_data()
    avEUb, clEUb, plEUb, tccEUb = load_euclideanBorder_data()
    integral_euclidean = compute_integral_euclidean(average_deg_quantile, small_world_quantile, avEU, clEU ./ plEU, length(velocities), length(thresholds))
    integral_euclidean_border = compute_integral_euclidean(average_deg_quantile, small_world_quantile, avEUb, clEUb ./ plEUb, length(velocities), length(thresholds))
    minimum_area, best_velocity, indices = find_minimum_v(integral_euclidean, velocities)
    println("Minimum area euclidean: ", minimum_area)
    println("Best velocity euclidean: ", best_velocity)
    println("Indices euclidean: ", indices)
    minimum_area, best_velocity, indices = find_minimum_v(integral_euclidean_border, velocities)
    println("Minimum area euclidean border: ", minimum_area)
    println("Best velocity euclidean border: ", best_velocity)
    println("Indices euclidean border: ", indices)
end
#CHANGE AND COMPUTE ONLY THINGS FOR SW

function myplot()
    subjects= readdlm("src/filtered-subjects-mod.txt", Int)
    Rrange = collect(0.5:0.5:14)
    thresholds = append!(collect(0.2:0.05:0.90), collect(0.92:0.02:0.98))
    average_deg_all, clustering_all, path_len_all, temp_clustering_all = load_fMRI_data(subjects, thresholds)
    average_deg_quantile, clustering_quantile, path_len_quantile, small_world_quantile, small_world_SB_quantile, temp_clustering_quantile = compute_quantiles(average_deg_all, clustering_all, path_len_all, temp_clustering_all, thresholds)
    average_deg_HY, clustering_HY, path_len_HY, temp_clustering_HY, small_world_HY, small_world_SB_HY = load_hyperbolic_data(length(Rrange))

    avEU, clEU, plEU, tccEU = load_euclidean_data()
    avEUb, clEUb, plEUb, tccEUb = load_euclideanBorder_data()
    p = plot(average_deg_quantile, small_world_SB_quantile[2, :], label="fMRI", xlabel="Average degree", ylabel="Small-worldness", legend=:topleft, xlims=(0, 170), ylims=(0, 0.6))
    plot!(p, average_deg_HY[8, :, 23], small_world_SB_HY[8, :, 23], label="HY")
    plot!(p, avEU[3,:], (tccEU ./ plEU)[3, :], label="EU")
    plot!(p, avEUb[4,:], (tccEUb ./ plEUb)[4, :], label="EU border")
end

main()