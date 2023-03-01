using NPZ, Random, Distributions,LinearAlgebra,Permutations,Plots,JLD2,Statistics, DelimitedFiles,Dierckx

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


alld=zeros(length(velocities),length(αrange))
for a in 1:length(αrange)
    @load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/hyperbolic_300_alpha_$(round(Int,αrange[a]*100))_R_05_05_14_v_01_01_09_correct.jld2" avHY clHY plHY tccHY 
    avHYct=zeros(length(velocities),length(thresholds))
    varvect=zeros(length(velocities),length(thresholds))
    avHYvect=avHY[:,:]
    varvect=clHY[:,:]./plHY[:,:]   #HERE
   
    for i in 1:length(velocities)
        spline=Spline1D(sort(reverse(avHYvect[i,:])),reverse(varvect[i,:]))
        splineorg=Spline1D(sort(reverse(avquatile[:])),reverse(smquantile[2,:]))
        range=collect(0:1:150)
        evspline=evaluate(spline,range)
        evsplineorg=evaluate(splineorg,range)
        diffspline=Spline1D(range,abs.(evspline-evsplineorg)) #HEREE
        alld[i,a]=integrate(diffspline,0,165)

    end    
    


end


println("Clustering smallworldness")
println(minimum(alld))
println(argmin(alld))
println(velocities[argmin(alld)[1]])
println(αrange[argmin(alld)[2]])

alld=zeros(length(velocities),length(αrange))
for a in 1:length(αrange)
    @load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/hyperbolic_300_alpha_$(round(Int,αrange[a]*100))_R_05_05_14_v_01_01_09_correct.jld2" avHY clHY plHY tccHY 
    avHYct=zeros(length(velocities),length(thresholds))
    varvect=zeros(length(velocities),length(thresholds))
    avHYvect=avHY[:,:]
    varvect=tccHY[:,:]./plHY[:,:]   #HERE
   
    for i in 1:length(velocities)
        spline=Spline1D(sort(reverse(avHYvect[i,:])),reverse(varvect[i,:]))
        splineorg=Spline1D(sort(reverse(avquatile[:])),reverse(smSBquantile[2,:]))
        range=collect(0:1:150)
        evspline=evaluate(spline,range)
        evsplineorg=evaluate(splineorg,range)
        diffspline=Spline1D(range,abs.(evspline-evsplineorg)) #HEREE
        alld[i,a]=integrate(diffspline,0,165)
    end    
    


end



p=plot(avquatile,smquantile[2,:],label="Real data 1050 Schaefer",lw=3)
for c in [(6,7),(9,20)]
    @load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/hyperbolic_300_alpha_$(round(Int,αrange[c[2]]*100))_R_05_05_14_v_01_01_09_correct.jld2" avHY clHY plHY tccHY 
    avHYvect=zeros(length(velocities),length(thresholds))
    varvect=zeros(length(velocities),length(thresholds))
    avHYvect=avHY[:,:]
    varvect=clHY[:,:]./plHY[:,:]
    p=plot!(avHYvect[c[1],:],varvect[c[1],:],label=" α=$(αrange[c[2]]), v=$(velocities[c[1]])",legendfontsize=6,legend=:topleft,lw=1,xlims=(0,165))
end



println("Sizemore and Basset smallworldness")
println(minimum(alld))
println(argmin(alld))
println(velocities[argmin(alld)[1]])
println(αrange[argmin(alld)[2]])