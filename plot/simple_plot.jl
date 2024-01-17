include("../optimization/findingParameters.jl")

function myplot()
    subjects = readdlm("src/filtered-subjects-mod.txt", Int)
    Rrange = collect(0.5:0.5:14)
    thresholds = append!(collect(0.2:0.05:0.90), collect(0.92:0.02:0.98))
    average_deg_all, clustering_all, path_len_all, temp_clustering_all = load_fMRI_data(subjects, thresholds)
    average_deg_quantile, clustering_quantile, path_len_quantile, small_world_quantile, small_world_SB_quantile, temp_clustering_quantile = compute_quantiles(average_deg_all, clustering_all, path_len_all, temp_clustering_all, thresholds)
    average_deg_HY, clustering_HY, path_len_HY, temp_clustering_HY, small_world_HY, small_world_SB_HY = load_hyperbolic_data(length(Rrange))

    avEU, clEU, plEU, tccEU = load_euclidean_data()
    avEUb, clEUb, plEUb, tccEUb = load_euclideanBorder_data()
    p = plot(average_deg_quantile, small_world_quantile[2, :], label="fMRI", xlabel="Average degree", ylabel="Small-worldness", legend=:topleft, xlims=(0, 170), ylims=(0, 0.6))
    plot!(p, average_deg_HY[1, :, 1], small_world_HY[1, :, 1], label="HY")
    plot!(p, avEU[1, :], (clEU./plEU)[1, :], label="EU")
    plot!(p, avEUb[9, :], (clEUb./plEUb)[9, :], label="EU border")
end

myplot()