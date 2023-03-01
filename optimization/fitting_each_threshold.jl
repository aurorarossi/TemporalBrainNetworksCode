using NPZ, Random, Distributions,LinearAlgebra,Permutations,Plots,JLD2,Statistics, DelimitedFiles,Dierckx, PlotlyJS

thresholds=append!(collect(0.2:0.05:0.90),collect(0.92:0.02:0.98))
subjects=readdlm("/user/aurossi/home/mri_networks/hyperbolicbrains/src/filtered-subjects-mod.txt",Int)
Rrange=collect(0.5:0.5:14)
radiusRange=collect(0.1:0.05:0.6)
avall=zeros(length(subjects),length(thresholds))
clall=zeros(length(subjects),length(thresholds))
plall=zeros(length(subjects),length(thresholds))
tccall=zeros(length(subjects),length(thresholds))


for i in 1:length(subjects)
    @load "/data/Hyperbrain/$(subjects[i])/$(subjects[i])_schaefer_300_thre_02_005_098_OR_RT_RE_ws60_wo30.jld2" av cl pl tcc avRT clRT plRT tccRT avRE clRE plRE tccRE
    avall[i,:]=av
    clall[i,:]=cl
    plall[i,:]=pl
    tccall[i,:]=tcc

end

avquatile=zeros(length(thresholds))
clquantile=zeros(3,length(thresholds))
plquantile=zeros(3,length(thresholds))
smquantile=zeros(3,length(thresholds))
smSBquantile=zeros(3,length(thresholds))
tccquantile=zeros(3,length(thresholds))

for j in 1:length(thresholds)
    avquatile[j]=mean(avall[:,j])
    clquantile[:,j]=quantile(clall[:,j])[2:4]
    plquantile[:,j]=quantile(plall[:,j])[2:4]
    smquantile[:,j]=quantile(clall[:,j]./plall[:,j])[2:4]
    smSBquantile[:,j]=quantile(tccall[:,j]./plall[:,j])[2:4]
    tccquantile[:,j]=quantile(tccall[:,j])[2:4]
    
end





αrange=collect(0.5:0.025:1.2)
velocities=collect(0.1:0.1:0.9)


alld=zeros(length(thresholds),length(velocities),length(αrange))
for a in 1:length(αrange)
    @load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/hyperbolic_300_alpha_$(round(Int,αrange[a]*100))_R_05_05_14_v_01_01_09_correct.jld2" avHY clHY plHY tccHY 
    avHYct=zeros(length(velocities),length(thresholds))
    varvect=zeros(length(velocities),length(thresholds))
    avHYvect=avHY[:,:]
    varvect=clHY[:,:]./plHY[:,:]   #HERE
   
    for i in 1:length(velocities)
        spline=Spline1D(sort(reverse(avHYvect[i,:])),reverse(varvect[i,:]))
        splineorg=Spline1D(sort(reverse(avquatile[:])),reverse(smquantile[2,:]))
        range=collect(0:1:165)
        evspline=evaluate(spline,range)
        evsplineorg=evaluate(splineorg,range)
        for j in 1:length(thresholds)
            alld[j,i,a]=abs.(evaluate(spline,avquatile[j])-evaluate(splineorg,avquatile[j]))
        end

    end    
    


end


alld_euclidean=zeros(length(velocities),length(thresholds))
@load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/euclidean_300_thre_02_005_098.jld2" avEU clEU plEU tccEU 
for i in 1:length(velocities)

        spline=Spline1D(sort(reverse(avEU[i,:])),reverse(clEU[i,:]./plEU[i,:]))
        splineorg=Spline1D(sort(reverse(avquatile[:])),reverse(smquantile[2,:]))
        range=collect(0:1:165)
        evspline=evaluate(spline,range)
        evsplineorg=evaluate(splineorg,range)
        for j in 1:length(thresholds)
            alld_euclidean[i,j]=abs.(evaluate(spline,avquatile[j])-evaluate(splineorg,avquatile[j]))
        end

end    
radiusRange=collect(0:0.05:0.8)

alld_euclideanBorder=zeros(length(velocities),length(radiusRange))
@load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/euclidean_300_thre_02_005_098.jld2" avEU clEU plEU tccEU 
for i in 1:length(velocities)

        spline=Spline1D(sort(avEU[i,:]),clEU[i,:]./plEU[i,:])
        splineorg=Spline1D(sort(reverse(avquatile[:])),reverse(smquantile[2,:]))
        range=collect(0:1:165)
        evspline=evaluate(spline,range)
        evsplineorg=evaluate(splineorg,range)
        for j in 1:length(radiusRange)
            alld_euclideanBorder[i,j]=abs.(evaluate(spline,avquatile[j])-evaluate(splineorg,avquatile[j]))
        end

end  
    

println("Clustering smallworldness")
for j in 1:length(thresholds)
    println(minimum(alld[j,:,:]))
    println(argmin(alld[j,:,:]))

end 


p=Plots.plot()
p=scatter!(avquatile,smquantile[2,:],label="Real data",color=:blue)
#p=scatter()

for j in 1:length(radiusRange)
    index_euB=argmin(alld_euclideanBorder[:,j])
 
    @load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/euclideanBorder_300_r_0_005_08.jld2" avEU clEU plEU tccEU 

    spline_euB=Spline1D(sort(avEU[index_euB,:]),clEU[index_euB,:]./plEU[index_euB,:])

    p=scatter!([avquatile[j]], [evaluate(spline_euB,avquatile[j])],color=:orange,alpha=.5,markershape=:hexagon,label="EuclideanBorder speed=$(velocities[index_euB])")


end 

for j in 1:length(thresholds)
    index=argmin(alld[j,:,:])
    index_eu=argmin(alld_euclidean[:,j])
    @load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/hyperbolic_300_alpha_$(round(Int,αrange[index[2]]*100))_R_05_05_14_v_01_01_09_correct.jld2" avHY clHY plHY tccHY 
    @load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/euclidean_300_thre_02_005_098.jld2" avEU clEU plEU tccEU 
    varvect=clHY[index[1],:]./plHY[index[1],:]
    spline=Spline1D(sort(reverse(avHY[index[1],:])),reverse(varvect))
    spline_eu=Spline1D(sort(reverse(avEU[index_eu,:])),reverse(clEU[index_eu,:]./plEU[index_eu,:]))
    p=scatter!([avquatile[j]], [evaluate(spline,avquatile[j])],color=:red,alpha=.5,markershape=:star5,label="Hyperbolic alpha=$(αrange[index[2]]),speed=$(velocities[index[1]])")
    p=scatter!([avquatile[j]], [evaluate(spline_eu,avquatile[j])],color=:green,alpha=.5,markershape=:rect,label="Euclidean speed=$(velocities[index_eu])",xlims=(0,165))
    

    
    xlabel!("Average degree")
    ylabel!("Small worldness")
    title!("Finding the set of parameter for each threshold")

end 
plotlyjs()
display(p)
plot!(size=(1300,1000))
Plots.png(p, "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/fittingeachthreshold")


pa=Plots.scatter()
for j in 1:length(thresholds)
    index=argmin(alld[j,:,:])
     pa=Plots.scatter!([avquatile[j]],[αrange[index[2]]],
    ylims=(0,1.4),
    xlims=(0,165),
    xlabel="Average degree",
    ylabel="Alpha",
    label="alpha=$(αrange[index[2]])")
end
display(pa)
Plots.html(pa, "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/bestalpha_allvelocities")


pv=Plots.scatter()
for j in 1:length(thresholds)
    index=argmin(alld[j,:,:])
    pv=Plots.scatter!([avquatile[j]],[velocities[index[1]]],
    ylims=(0,1),
    xlims=(0,165),
    xlabel="Thresholds",
    ylabel="Velocities",
    label="velocities=$(velocities[index[1]])")
end
display(pv)
Plots.html(pv, "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/bestvelocities_allalpha")


