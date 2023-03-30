using NPZ, Random, Distributions, LinearAlgebra, Permutations, Plots, JLD2, Statistics, DelimitedFiles, Interpolations, QuadGK, LaTeXStrings, Measures

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
        @load "/user/aurossi/home/mri_networks/TemporalBrainNetworksCode/data/euclidean_300_thre_02_005_098_v01_$(i).jld2" avEU clEU plEU tccEU
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
        @load "/user/aurossi/home/mri_networks/TemporalBrainNetworksCode/data/euclideanBorder_300_r_0_005_08_v01_$(i).jld2" avEU clEU plEU tccEU
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
        @load "/user/aurossi/home/mri_networks/TemporalBrainNetworksCode/data/hyperbolic_300_alpha_50_R_0_05_18_v_10_$(i).jld2" avHY clHY plHY tccHY
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
        @load "/user/aurossi/home/mri_networks/TemporalBrainNetworksCode/data/hyperbolic_300_alpha_105_R_0_05_18_v_80_$(i).jld2" avHY clHY plHY tccHY
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
#L"\mathrm{Temporal\,\, Small\,\, Worldness}"
function myplot_withmarker(measure, variable, variableRT, variableRE, variableEU, variableEUB, variableHY, variableHYSB, av, avRT, avRE, avEU, avEUB, avHY, avHYSB)
    p = plot(av, variable[1, :], fillrange=variable[3, :], fillalpha=0.35, c=1, label="", lw=0, palette=:seaborn_colorblind)
    plot!(av, variable[2, :], label=L"\mathrm{Real\,\, data\,\, subjects\,\, Schaefer}", markershape=:circle, lw=2, c=1)
    plot!(avRE, variableRE[1, :], fillrange=variableRE[3, :], fillalpha=0.35, c=2, label="", lw=0)
    plot!(avRE, variableRE[2, :], label=L"\mathrm{Random \,\, Edges}", markershape=:pentagon, lw=2, c=2)
    plot!(avRT, variableRT[1, :], fillrange=variableRT[3, :], fillalpha=0.35, c=3, label="", lw=0)
    plot!(avRT, variableRT[2, :], label=L"\mathrm{Random \,\,Times}", markershape=:hexagon, lw=2, c=3)
    plot!(avEU, variableEU[1, :], fillrange=variableEU[3, :], fillalpha=0.35, c=4, label="", lw=0)
    plot!(avEU, variableEU[2, :], label=L"\mathrm{Geometric \,\,on\,\, torus\,\, }v=0.1", markershape=:rect, lw=2, c=4)
    plot!(avEUB, variableEUB[1, :], fillrange=variableEUB[3, :], fillalpha=0.35, c=5, label="", lw=0)
    plot!(avEUB, variableEUB[2, :], label=L"\mathrm{Geometric \,\,}v=0.1", markershape=:rect, lw=2, c=5, legend=:topleft, legendfontsize=5, xlims=(0, 165), grid=false, dpi=1200)

    plot!(avHY, variableHY[1, :], fillrange=variableHY[3, :], fillalpha=0.35, c=10, label="", lw=0)
    plot!(avHY, variableHY[2, :], label=L"\mathrm{Hyperbolic \,\, }α=0.5, v=0.1, ζ=1", c=10, lw=1, markershape=:star8)
    plot!(avHYSB, variableHYSB[1, :], fillrange=variableHYSB[3, :], fillalpha=0.35, c=9, label="", lw=0)
    plot!(avHYSB, variableHYSB[2, :], label=L"\mathrm{Hyperbolic \,\, }α=1.05, v=0.8, ζ=1 ", c=9, lw=1, markershape=:star8)

    xlabel!(L"\mathrm{Average\,\, degree}")
    ylabel!("$(measure)")
    title!(L"\mathrm{302\,\, nodes}")
    return p
end

function myplot(measure, legendplace, yaxisscale, xaxisscale, variable, variableRT, variableRE, variableEU, variableEUB, variableHY, variableHYSB, av, avRT, avRE, avEU, avEUB, avHY, avHYSB)
    line_lw = 3
    fillalpha = 0.35
    red = 1
    blue = 2
    green = 3
    pink = 8
    orange = 5
    violet = 4
    grey = 9
    if measure !="Temporal Path Length"
        yticks = ([0, 0.25, 0.5,0.75,1], [L"0", L"0.25", L"0.5", L"0.75",L"1"])
    else
        yticks = ([10^0,10^1], [L"10^0", L"10^1"])
    end
    p = plot(av, variable[1, :], fillrange=variable[3, :], fillalpha=fillalpha, c=blue, label="", lw=0, palette=:Set1_9, legend=legendplace, legendfontsize=12, xlims=(0, 170), grid=false, dpi=1200, yaxis=yaxisscale, xaxis=xaxisscale, xticks = ([0, 50, 100, 150, 200], [L"0", L"50", L"100", L"150", L"200"]), yticks = yticks)

    plot!(avRE, variableRE[1, :], fillrange=variableRE[3, :], fillalpha=fillalpha, c=pink, label="", lw=0)
    plot!(avRE, variableRE[2, :], label=L"\mathrm{Random \,\, Edges}", lw=line_lw, c=pink)
    plot!(avRT, variableRT[1, :], fillrange=variableRT[3, :], fillalpha=fillalpha, c=grey, label="", lw=0)
    plot!(avRT, variableRT[2, :], label=L"\mathrm{Random \,\,Times}", lw=line_lw, c=grey)
    plot!(avEU, variableEU[1, :], fillrange=variableEU[3, :], fillalpha=fillalpha, c=orange, label="", lw=0)
    plot!(avEU, variableEU[2, :], label=L"\mathrm{Geometric \,\,on\,\, torus\,\, }v=0.1", lw=line_lw, c=orange)
    plot!(avEUB, variableEUB[1, :], fillrange=variableEUB[3, :], fillalpha=fillalpha, c=green, label="", lw=0)
    plot!(avEUB, variableEUB[2, :], label=L"\mathrm{Geometric \,\,}v=0.1", lw=line_lw, c=green,)
    plot!(avHY, variableHY[1, :], fillrange=variableHY[3, :], fillalpha=fillalpha, c=violet, label="", lw=0)
    plot!(avHY, variableHY[2, :], label=L"\mathrm{Hyperbolic \,\, }α=0.5, v=0.1, ζ=1", c=violet, lw=line_lw)
    plot!(avHYSB, variableHYSB[1, :], fillrange=variableHYSB[3, :], fillalpha=fillalpha, c=red, label="", lw=0)
    plot!(avHYSB, variableHYSB[2, :], label=L"\mathrm{Hyperbolic \,\, }α=1.05, v=0.8, ζ=1 ", c=red, lw=line_lw)
    plot!(av, variable[2, :], label=L"\mathrm{Real\,\, data\,\,  Schaefer}", lw=line_lw, c=blue)

    xlabel!(L"\mathrm{Average\,\,  degree}")
    if measure == "Temporal Clustering"
        ylabel!(L"\mathrm{Temporal\,\,clustering \,\, coefficient \,\,} C")
    elseif measure == "Temporal Path Length"
        ylabel!(L"\log_{10}(\mathrm{Temporal \,\, Path \,\, length})\,\, \log_{10}(L)")
    elseif measure == "Temporal Small Worldness"
        ylabel!(L"\mathrm{Temporal\,\, small\,\, worldness \,\,} S")
    elseif measure == "Temporal Small Worldness S_{SB}"
        ylabel!(L"\mathrm{Temporal\,\, small\,\, worldness \,\,} S_{SB}")
    elseif measure == "Temporal Correlation Coefficient"
        ylabel!(L"\mathrm{Temporal\,\, correlation \,\, coefficient \,\,} TC")
    end

    return p
end

function main()
    average_deg_quantile, clustering_quantile, path_len_quantile, small_world_quantile, small_world_SB_quantile, temp_corr_quantile, average_deg_quantile_RT, clustering_quantile_RT, path_len_quantile_RT, small_world_quantile_RT, small_world_SB_quantile_RT, temp_corr_quantile_RT, average_deg_quantile_RE, clustering_quantile_RE, path_len_quantile_RE, small_world_quantile_RE, small_world_SB_quantile_RE, temp_corr_quantile_RE, average_deg_quantile_EU, clustering_quantile_EU, path_len_quantile_EU, small_world_quantile_EU, small_world_SB_quantile_EU, temp_corr_quantile_EU, average_deg_quantile_EUB, clustering_quantile_EUB, path_len_quantile_EUB, small_world_quantile_EUB, small_world_SB_quantile_EUB, temp_corr_quantile_EUB, average_deg_quantile_HY, clustering_quantile_HY, path_len_quantile_HY, small_world_quantile_HY, small_world_SB_quantile_HY, temp_corr_quantile_HY, average_deg_quantile_HYSB, clustering_quantile_HYSB, path_len_quantile_HYSB, small_world_quantile_HYSB, small_world_SB_quantile_HYSB, temp_corr_quantile_HYSB = load_and_quantiles()
    p1 = myplot("Temporal Small Worldness", :topleft, :identity, :identity, small_world_quantile, small_world_quantile_RT, small_world_quantile_RE, small_world_quantile_EU, small_world_quantile_EUB, small_world_quantile_HY, small_world_quantile_HYSB, average_deg_quantile, average_deg_quantile_RT, average_deg_quantile_RE, average_deg_quantile_EU, average_deg_quantile_EUB, average_deg_quantile_HY, average_deg_quantile_HYSB)
    p2 = myplot("Temporal Small Worldness S_{SB}", false, :identity, :identity, small_world_SB_quantile, small_world_SB_quantile_RT, small_world_SB_quantile_RE, small_world_SB_quantile_EU, small_world_SB_quantile_EUB, small_world_SB_quantile_HY, small_world_SB_quantile_HYSB, average_deg_quantile, average_deg_quantile_RT, average_deg_quantile_RE, average_deg_quantile_EU, average_deg_quantile_EUB, average_deg_quantile_HY, average_deg_quantile_HYSB)
    p3 = myplot("Temporal Correlation Coefficient", false, :identity, :identity, temp_corr_quantile, temp_corr_quantile_RT, temp_corr_quantile_RE, temp_corr_quantile_EU, temp_corr_quantile_EUB, temp_corr_quantile_HY, temp_corr_quantile_HYSB, average_deg_quantile, average_deg_quantile_RT, average_deg_quantile_RE, average_deg_quantile_EU, average_deg_quantile_EUB, average_deg_quantile_HY, average_deg_quantile_HYSB)
    p4 = myplot("Temporal Clustering", false, :identity, :identity, clustering_quantile, clustering_quantile_RT, clustering_quantile_RE, clustering_quantile_EU, clustering_quantile_EUB, clustering_quantile_HY, clustering_quantile_HYSB, average_deg_quantile, average_deg_quantile_RT, average_deg_quantile_RE, average_deg_quantile_EU, average_deg_quantile_EUB, average_deg_quantile_HY, average_deg_quantile_HYSB)
    p5 = myplot("Temporal Path Length", false, :log, :identity, path_len_quantile, path_len_quantile_RT, path_len_quantile_RE, path_len_quantile_EU, path_len_quantile_EUB, path_len_quantile_HY, path_len_quantile_HYSB, average_deg_quantile, average_deg_quantile_RT, average_deg_quantile_RE, average_deg_quantile_EU, average_deg_quantile_EUB, average_deg_quantile_HY, average_deg_quantile_HYSB)
    p = plot(p1, p2, p3, p4, p5, layout=(5, 1), size=(800, 3000), margin=25mm)
    display(p1)
    display(p2)
    display(p3)
    display(p4)
    display(p5)
    savefig(p1, "/user/aurossi/home/mri_networks/TemporalBrainNetworksCode/images/tsw.png")
    savefig(p2, "/user/aurossi/home/mri_networks/TemporalBrainNetworksCode/images/tswsb.png")
    savefig(p3, "/user/aurossi/home/mri_networks/TemporalBrainNetworksCode/images/tcc.png")
    savefig(p4, "/user/aurossi/home/mri_networks/TemporalBrainNetworksCode/images/tclustering.png")
    savefig(p5, "/user/aurossi/home/mri_networks/TemporalBrainNetworksCode/images/tpath.png")
end

main()