using NPZ, Random, Distributions, LinearAlgebra, Permutations, Plots, JLD2, Statistics, DelimitedFiles, Interpolations, QuadGK, LaTeXStrings, Measures, Dierckx

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






function plot_hyperbolic(thresholds,velocities,αrange, average_deg_quantile, small_world_quantile, p)
    alld=zeros(length(thresholds),length(velocities),length(αrange))
    for a in 1:length(αrange)
        @load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/hyperbolic_300_alpha_$(round(Int,αrange[a]*100))_R_05_05_14_v_01_01_09_correct.jld2" avHY clHY plHY tccHY 
        avHYvect=zeros(length(velocities),length(thresholds))
        varvect=zeros(length(velocities),length(thresholds))
        avHYvect=avHY[:,:]
        varvect=clHY[:,:]./plHY[:,:]   #HERE
        for i in 1:length(velocities)
            spline=Spline1D(sort(reverse(avHYvect[i,:])),reverse(varvect[i,:]))
            splineorg=Spline1D(sort(reverse(average_deg_quantile[:])),reverse(small_world_quantile[2,:]))
            for j in 1:length(thresholds)
                alld[j,i,a]=abs.(evaluate(spline,average_deg_quantile[j])-evaluate(splineorg,average_deg_quantile[j]))
            end

        end    
    end
    for j in 1:length(thresholds)-1
        index=argmin(alld[j,:,:])
        @load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/hyperbolic_300_alpha_$(round(Int,αrange[index[2]]*100))_R_05_05_14_v_01_01_09_correct.jld2" avHY clHY plHY tccHY 
        varvect=clHY[index[1],:]./plHY[index[1],:]
        spline=Spline1D(sort(reverse(avHY[index[1],:])),reverse(varvect))
        scatter!(p,[average_deg_quantile[j]], [evaluate(spline,average_deg_quantile[j])],color=1,alpha=0.75,markersize=4,markerstrokewidth=0.6,markershape=:diamond,label="")
    end
    index=argmin(alld[end,:,:])
    @load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/hyperbolic_300_alpha_$(round(Int,αrange[index[2]]*100))_R_05_05_14_v_01_01_09_correct.jld2" avHY clHY plHY tccHY 
    varvect=clHY[index[1],:]./plHY[index[1],:]
    spline=Spline1D(sort(reverse(avHY[index[1],:])),reverse(varvect))
    scatter!(p,[average_deg_quantile[end]], [evaluate(spline,average_deg_quantile[end])],color=1,alpha=0.75,markersize=4,markerstrokewidth=0.6,markershape=:diamond, label=L"\mathrm{RTH}")

end

function plot_euclidean(thresholds,velocities, average_deg_quantile, small_world_quantile,p)
    alld_euclidean=zeros(length(velocities),length(thresholds))
    @load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/euclidean_300_thre_02_005_098.jld2" avEU clEU plEU tccEU 
    for i in 1:length(velocities)
            spline=Spline1D(sort(reverse(avEU[i,:])),reverse(clEU[i,:]./plEU[i,:]))
            splineorg=Spline1D(sort(reverse(average_deg_quantile[:])),reverse(small_world_quantile[2,:]))
            for j in 1:length(thresholds)
                alld_euclidean[i,j]=abs.(evaluate(spline,average_deg_quantile[j])-evaluate(splineorg,average_deg_quantile[j]))
            end
    end    
    for j in 1:length(thresholds)-1
        index_eu=argmin(alld_euclidean[:,j])
        @load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/euclidean_300_thre_02_005_098.jld2" avEU clEU plEU tccEU 
        spline_eu=Spline1D(sort(reverse(avEU[index_eu,:])),reverse(clEU[index_eu,:]./plEU[index_eu,:]))
        scatter!(p,[average_deg_quantile[j]], [evaluate(spline_eu,average_deg_quantile[j])],color=5,alpha=0.75,markersize=4,markerstrokewidth=0.6,markershape=:rect,label="")
    end
    index_eu=argmin(alld_euclidean[:,end])
    @load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/euclidean_300_thre_02_005_098.jld2" avEU clEU plEU tccEU 
    spline_eu=Spline1D(sort(reverse(avEU[index_eu,:])),reverse(clEU[index_eu,:]./plEU[index_eu,:]))
    scatter!(p,[average_deg_quantile[end]], [evaluate(spline_eu,average_deg_quantile[end])],color=5,alpha=0.75,markersize=4,markerstrokewidth=0.6,markershape=:rect, label= L"\mathrm{RTT}")
end

function plot_euclideanBorder(thresholds,velocities,radiusRange, average_deg_quantile, small_world_quantile,p)
    alld_euclideanBorder=zeros(length(velocities),length(radiusRange))
    @load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/euclideanBorder_300_r_0_005_08.jld2" avEU clEU plEU tccEU 
    for i in 1:length(velocities)
            spline=Spline1D(sort(avEU[i,:]),clEU[i,:]./plEU[i,:])
            splineorg=Spline1D(sort(reverse(average_deg_quantile[:])),reverse(small_world_quantile[2,:]))
            for j in 1:length(radiusRange)
                alld_euclideanBorder[i,j]=abs.(evaluate(spline,average_deg_quantile[j])-evaluate(splineorg,average_deg_quantile[j]))
            end
    end  
    for j in 1:length(radiusRange)-1
        index_euB=argmin(alld_euclideanBorder[:,j])
        @load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/euclideanBorder_300_r_0_005_08.jld2" avEU clEU plEU tccEU 
        spline_euB=Spline1D(sort(avEU[index_euB,:]),clEU[index_euB,:]./plEU[index_euB,:])
        scatter!(p,[average_deg_quantile[j]], [evaluate(spline_euB,average_deg_quantile[j])],color=3,alpha=0.75,markersize=4,markerstrokewidth=0.6,markershape=:hexagon,label="")
    end 
    index_euB=argmin(alld_euclideanBorder[:,end])
    @load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/euclideanBorder_300_r_0_005_08.jld2" avEU clEU plEU tccEU 
    spline_euB=Spline1D(sort(avEU[index_euB,:]),clEU[index_euB,:]./plEU[index_euB,:])
    scatter!(p,[average_deg_quantile[end]], [evaluate(spline_euB,average_deg_quantile[end])],color=3,alpha=0.75,markersize=4,markerstrokewidth=0.6,markershape=:hexagon,label=L"\mathrm{RTS}")
end


function main()
    velocities=collect(0.1:0.1:0.9)
    thresholds=append!(collect(0.2:0.05:0.90),collect(0.92:0.02:0.98))
    radiusRange=collect(0:0.05:0.8)
    αrange=collect(0.5:0.025:1.2)

    subjects = readdlm("src/filtered-subjects-mod.txt", Int)
    average_deg_all, clustering_all, path_len_all, temp_corr_all = load_fMRI_data(subjects, thresholds)
    average_deg_quantile, clustering_quantile, path_len_quantile, small_world_quantile, small_world_SB_quantile, temp_corr_quantile = compute_quantiles(average_deg_all, clustering_all, path_len_all, temp_corr_all, length(thresholds))

    p = Plots.scatter(average_deg_quantile,small_world_quantile[2,:],c=2, palette=:Set1_9, legendfontsize=12, xlims=(0, 170), grid=false,xticks = ([0, 50, 100, 150, 200], [L"0", L"50", L"100", L"150", L"200"]), yticks = ([ 0.25, 0.5,0.75,1], [ L"0.25", L"0.5", L"0.75",L"1"]),markersize=5.3,markerstrokewidth=0.6,label=L"\mathrm{Real\,\, data\,\,  Schaefer}", dpi=1200)
    ylabel!(L"\mathrm{Temporal\,\, small\,\, worldness \,\,} S")
    xlabel!(L"\mathrm{Average\,\,  degree}")
    plot_euclideanBorder(thresholds,velocities,radiusRange, average_deg_quantile, small_world_quantile,p)
    plot_hyperbolic(thresholds,velocities,αrange, average_deg_quantile, small_world_quantile, p)
    plot_euclidean(thresholds,velocities, average_deg_quantile, small_world_quantile,p)
    display(p)
    savefig(p, "/user/aurossi/home/mri_networks/TemporalBrainNetworksCode/images/swEachThreshold.png")
end

main()