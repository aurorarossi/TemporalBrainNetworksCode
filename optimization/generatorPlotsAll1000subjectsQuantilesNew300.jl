using NPZ, Random, Distributions,LinearAlgebra,Permutations,Plots,JLD2,Statistics, DelimitedFiles, LaTeXStrings


 

thresholds=append!(collect(0.2:0.05:0.90),collect(0.92:0.02:0.98))
#CREATION
subjects=readdlm("/user/aurossi/home/mri_networks/hyperbolicbrains/src/filtered-subjects-mod.txt",Int)
Rrange=collect(0.5:0.5:14)
radiusRange=collect(0:0.05:0.6)
avall=zeros(length(subjects),length(thresholds))
clall=zeros(length(subjects),length(thresholds))
plall=zeros(length(subjects),length(thresholds))
tccall=zeros(length(subjects),length(thresholds))
avallRT=zeros(length(subjects),length(thresholds))
clallRT=zeros(length(subjects),length(thresholds))
plallRT=zeros(length(subjects),length(thresholds))
tccallRT=zeros(length(subjects),length(thresholds))
avallRE=zeros(length(subjects),length(thresholds))
clallRE=zeros(length(subjects),length(thresholds))
plallRE=zeros(length(subjects),length(thresholds))
tccallRE=zeros(length(subjects),length(thresholds))
avallEU=zeros(10,length(thresholds))
clallEU=zeros(10,length(thresholds))
plallEU=zeros(10,length(thresholds))
tccallEU=zeros(10,length(thresholds))
avallEUB=zeros(10,length(radiusRange))
clallEUB=zeros(10,length(radiusRange))
plallEUB=zeros(10,length(radiusRange))
tccallEUB=zeros(10,length(radiusRange))

avallHYsm=zeros(10,length(Rrange))
clallHYsm=zeros(10,length(Rrange))
plallHYsm=zeros(10,length(Rrange))
tccallHYsm=zeros(10,length(Rrange))
avallHYsmSB=zeros(10,length(Rrange))
clallHYsmSB=zeros(10,length(Rrange))
plallHYsmSB=zeros(10,length(Rrange))
tccallHYsmSB=zeros(10,length(Rrange))

for i in 1:length(subjects)
    @load "/data/Hyperbrain/$(subjects[i])/$(subjects[i])_schaefer_300_thre_02_005_098_OR_RT_RE_ws60_wo30.jld2" av cl pl tcc
    @load "/data/Hyperbrain/$(subjects[i])/$(subjects[i])_schaefer_300_thre_02_005_098_permuted_ws60_wo30.jld2" avRT clRT plRT tccRT avRE clRE plRE tccRE
    avall[i,:]=av
    clall[i,:]=cl
    plall[i,:]=pl
    tccall[i,:]=tcc
    avallRT[i,:]=avRT
    clallRT[i,:]=clRT
    plallRT[i,:]=plRT
    tccallRT[i,:]=tccRT
    avallRE[i,:]=avRE
    clallRE[i,:]=clRE
    plallRE[i,:]=plRE
    tccallRE[i,:]=tccRE
end
for i in 1:10
    # @load "/user/aurossi/home/mri_networks/hyperbolic_300_alpha_88_R_05_05_10_v_01_$(i).jld2" avHY clHY plHY tccHY
    # avallHY[i,:]=avHY[1,:]
    # clallHY[i,:]=clHY[1,:]
    # plallHY[i,:]=plHY[1,:]
    # tccallHY[i,:]=tccHY[1,:]
    @load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/hyperbolic_300_alpha_65_R_05_05_14_v_60_$(i)_sm.jld2" avHY clHY plHY tccHY
    avallHYsm[i,:]=avHY[1,:]
    clallHYsm[i,:]=clHY[1,:]
    plallHYsm[i,:]=plHY[1,:]
    tccallHYsm[i,:]=tccHY[1,:]
    @load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/hyperbolic_300_alpha_98_R_05_05_14_v_90_$(i)_sm.jld2" avHY clHY plHY tccHY
    avallHYsmSB[i,:]=avHY[1,:]
    clallHYsmSB[i,:]=clHY[1,:]
    plallHYsmSB[i,:]=plHY[1,:]
    tccallHYsmSB[i,:]=tccHY[1,:]
    @load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/euclidean_300_thre_02_005_098_v_09_$(i).jld2" avEU clEU plEU tccEU 
    avallEU[i,:]=avEU
    clallEU[i,:]=clEU
    plallEU[i,:]=plEU
    tccallEU[i,:]=tccEU

    @load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/euclideanBorder_300_r_01_005_08_v_09_$(i).jld2" avEU clEU plEU tccEU
    avallEUB[i,:]=avEU
    clallEUB[i,:]=clEU
    plallEUB[i,:]=plEU
    tccallEUB[i,:]=tccEU
end
avquatile=zeros(length(thresholds))
clquantile=zeros(3,length(thresholds))
plquantile=zeros(3,length(thresholds))
smquantile=zeros(3,length(thresholds))
smSBquantile=zeros(3,length(thresholds))
tccquantile=zeros(3,length(thresholds))
avquatileRT=zeros(length(thresholds))
clquantileRT=zeros(3,length(thresholds))
plquantileRT=zeros(3,length(thresholds))
smquantileRT=zeros(3,length(thresholds))
smSBquantileRT=zeros(3,length(thresholds))
tccquantileRT=zeros(3,length(thresholds))
avquatileRE=zeros(length(thresholds))
clquantileRE=zeros(3,length(thresholds))
plquantileRE=zeros(3,length(thresholds))
smquantileRE=zeros(3,length(thresholds))
smSBquantileRE=zeros(3,length(thresholds))
tccquantileRE=zeros(3,length(thresholds))
avquatileEU=zeros(length(thresholds))
clquantileEU=zeros(3,length(thresholds))
plquantileEU=zeros(3,length(thresholds))
smquantileEU=zeros(3,length(thresholds))
smSBquantileEU=zeros(3,length(thresholds))
tccquantileEU=zeros(3,length(thresholds))
avquatileEUB=zeros(length(radiusRange))
clquantileEUB=zeros(3,length(radiusRange))
plquantileEUB=zeros(3,length(radiusRange))
smquantileEUB=zeros(3,length(radiusRange))
smSBquantileEUB=zeros(3,length(radiusRange))
tccquantileEUB=zeros(3,length(radiusRange))

avquatileHYsm=zeros(length(Rrange))
clquantileHYsm=zeros(3,length(Rrange))
plquantileHYsm=zeros(3,length(Rrange))
smquantileHYsm=zeros(3,length(Rrange))
smSBquantileHYsm=zeros(3,length(Rrange))
tccquantileHYsm=zeros(3,length(Rrange))

avquatileHYsmSB=zeros(length(Rrange))
clquantileHYsmSB=zeros(3,length(Rrange))
plquantileHYsmSB=zeros(3,length(Rrange))
smquantileHYsmSB=zeros(3,length(Rrange))
smSBquantileHYsmSB=zeros(3,length(Rrange))
tccquantileHYsmSB=zeros(3,length(Rrange))

for j in 1:length(thresholds)
    avquatile[j]=mean(avall[:,j])
    clquantile[:,j]=quantile(clall[:,j])[2:4]
    plquantile[:,j]=quantile(plall[:,j])[2:4]
    smquantile[:,j]=quantile(clall[:,j]./plall[:,j])[2:4]
    smSBquantile[:,j]=quantile(tccall[:,j]./plall[:,j])[2:4]
    tccquantile[:,j]=quantile(tccall[:,j])[2:4]
    avquatileRT[j]=mean(avallRT[:,j])
    clquantileRT[:,j]=quantile(clallRT[:,j])[2:4]
    plquantileRT[:,j]=quantile(plallRT[:,j])[2:4]
    smquantileRT[:,j]=quantile(clallRT[:,j]./plallRT[:,j])[2:4]
    smSBquantileRT[:,j]=quantile(tccallRT[:,j]./plallRT[:,j])[2:4]
    tccquantileRT[:,j]=quantile(tccallRT[:,j])[2:4]
    avquatileRE[j]=mean(avallRE[:,j])
    clquantileRE[:,j]=quantile(clallRE[:,j])[2:4]
    plquantileRE[:,j]=quantile(plallRE[:,j])[2:4]
    smquantileRE[:,j]=quantile(clallRE[:,j]./plallRE[:,j])[2:4]
    smSBquantileRE[:,j]=quantile(tccallRE[:,j]./plallRE[:,j])[2:4]
    tccquantileRE[:,j]=quantile(tccallRE[:,j])[2:4]
    avquatileEU[j]=mean(avallEU[:,j])
    clquantileEU[:,j]=quantile(clallEU[:,j])[2:4]
    plquantileEU[:,j]=quantile(plallEU[:,j])[2:4]
    smquantileEU[:,j]=quantile(clallEU[:,j]./plallEU[:,j])[2:4]
    smSBquantileEU[:,j]=quantile(tccallEU[:,j]./plallEU[:,j])[2:4]
    tccquantileEU[:,j]=quantile(tccallEU[:,j])[2:4]
end
for j in 1:length(radiusRange)
    avquatileEUB[j]=mean(avallEUB[:,j])
    clquantileEUB[:,j]=quantile(clallEUB[:,j])[2:4]
    plquantileEUB[:,j]=quantile(plallEUB[:,j])[2:4]
    smquantileEUB[:,j]=quantile(clallEUB[:,j]./plallEUB[:,j])[2:4]
    smSBquantileEUB[:,j]=quantile(tccallEUB[:,j]./plallEUB[:,j])[2:4]
    tccquantileEUB[:,j]=quantile(tccallEUB[:,j])[2:4]
end
for j in 1:length(Rrange)
    avquatileHYsm[j]=mean(avallHYsm[:,j])
    clquantileHYsm[:,j]=quantile(clallHYsm[:,j])[2:4]
    plquantileHYsm[:,j]=quantile(plallHYsm[:,j])[2:4]
    smquantileHYsm[:,j]=quantile(clallHYsm[:,j]./plallHYsm[:,j])[2:4]
    smSBquantileHYsm[:,j]=quantile(tccallHYsm[:,j]./plallHYsm[:,j])[2:4]
    tccquantileHYsm[:,j]=quantile(tccallHYsm[:,j])[2:4]

    avquatileHYsmSB[j]=mean(avallHYsmSB[:,j])
    clquantileHYsmSB[:,j]=quantile(clallHYsmSB[:,j])[2:4]
    plquantileHYsmSB[:,j]=quantile(plallHYsmSB[:,j])[2:4]
    smquantileHYsmSB[:,j]=quantile(clallHYsmSB[:,j]./plallHYsmSB[:,j])[2:4]
    smSBquantileHYsmSB[:,j]=quantile(tccallHYsmSB[:,j]./plallHYsmSB[:,j])[2:4]
    tccquantileHYsmSB[:,j]=quantile(tccallHYsmSB[:,j])[2:4]
    
end

#SMALL WORLDNESS
    
p=plot(avquatile,smquantile[1,:],fillrange=smquantile[3,:],fillalpha = 0.35, c = 1, label="",lw = 0,palette=:seaborn_colorblind)
plot!(avquatile,smquantile[2,:],label=L"\mathrm{Real\,\, data\,\, subjects\,\, Schaefer}",markershape=:circle,lw = 2,c=1)
plot!(avquatile,smquantileRE[1,:],fillrange=smquantileRE[3,:],fillalpha = 0.35, c=2, label="",lw = 0)
plot!(avquatile,smquantileRE[2,:],label=L"\mathrm{Random \,\, Edges}",markershape=:pentagon,lw = 2,c=2)
plot!(avquatile,smquantileRT[1,:],fillrange=smquantileRT[3,:],fillalpha = 0.35, c=3, label="",lw = 0)
plot!(avquatile,smquantileRT[2,:],label=L"\mathrm{Random \,\,Times}",markershape=:hexagon,lw = 2,c=3)
plot!(avquatileEU,smquantileEU[1,:],fillrange=smquantileEU[3,:],fillalpha = 0.35, c=4, label="",lw = 0)
plot!(avquatileEU,smquantileEU[2,:],label=L"\mathrm{Geometric \,\,on\,\, torus\,\, }v=0.9",markershape=:rect,lw = 2,c=4)
plot!(avquatileEUB,smquantileEUB[1,:],fillrange=smquantileEUB[3,:],fillalpha = 0.35, c=5, label="",lw = 0)
plot!(avquatileEUB,smquantileEUB[2,:],label=L"\mathrm{Geometric \,\,}v=0.9",markershape=:rect,lw = 2,c=5,legend=:topleft,legendfontsize=5,xlims=(0, 165),grid=false,dpi=1200)

plot!(avquatileHYsm,smquantileHYsm[1,:],fillrange=smquantileHYsm[3,:],fillalpha = 0.35, c=10, label="",lw = 0)
plot!(avquatileHYsm,smquantileHYsm[2,:],label=L"\mathrm{Hyperbolic \,\, }α=1.120, v=0.6, ζ=1 \mathrm{\,\, (SW)}",c=10,lw = 1,markershape=:star8)
plot!(avquatileHYsmSB,smquantileHYsmSB[1,:],fillrange=smquantileHYsmSB[3,:],fillalpha = 0.35, c=9, label="",lw = 0)
plot!(avquatileHYsmSB,smquantileHYsmSB[2,:],label=L"\mathrm{Hyperbolic \,\, }α=1.125, v=0.5, ζ=1 \mathrm{\,\, (SWSB)}",c=9,lw = 1,markershape=:star8)

    xlabel!(L"\mathrm{Average\,\, degree}")
    ylabel!(L"\mathrm{Temporal\,\, Small\,\, Worldness}")
    title!(L"\mathrm{302\,\, nodes}")
    png("/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/optimizatedSmallWorldness300")

# #CLUSTERING

p1=plot(avquatile,clquantile[1,:],fillrange=clquantile[3,:],fillalpha = 0.35, c = 1, label="",lw = 0,palette=:seaborn_colorblind)
plot!(avquatile,clquantile[2,:],label=L"\mathrm{Real\,\, data\,\, subjects\,\, Schaefer}",markershape=:circle,lw = 2,c=1)
plot!(avquatile,clquantileRE[1,:],fillrange=clquantileRE[3,:],fillalpha = 0.35, c=2, label="",lw = 0)
plot!(avquatile,clquantileRE[2,:],label=L"\mathrm{Random \,\, Edges}",markershape=:pentagon,lw = 2,c=2)
plot!(avquatile,clquantileRT[1,:],fillrange=clquantileRT[3,:],fillalpha = 0.35, c=3, label="",lw = 0)
plot!(avquatile,clquantileRT[2,:],label=L"\mathrm{Random \,\,Times}",markershape=:hexagon,lw = 2,c=3)
plot!(avquatileEU,clquantileEU[1,:],fillrange=clquantileEU[3,:],fillalpha = 0.35, c=4, label="",lw = 0)
plot!(avquatileEU,clquantileEU[2,:],label=L"\mathrm{Geometric \,\,on\,\, torus\,\, v=0.9}",markershape=:rect,lw = 2,c=4)
plot!(avquatileEUB,clquantileEUB[1,:],fillrange=clquantileEUB[3,:],fillalpha = 0.35, c=5, label="",lw = 0)
plot!(avquatileEUB,clquantileEUB[2,:],label=L"\mathrm{Geometric\,\,} v=0.9",markershape=:rect,lw = 2,c=5,legend=:bottomright,legendfontsize=5,xlims=(0, 165),grid=false,dpi=1200)
plot!(avquatileHYsm,clquantileHYsm[1,:],fillrange=clquantileHYsm[3,:],fillalpha = 0.35, c=10, label="",lw = 0)
plot!(avquatileHYsm,clquantileHYsm[2,:],label=L"\mathrm{Hyperbolic \,\, }α=1.120, v=0.6, ζ=1 \mathrm{\,\, (SW)}",c=10,lw = 1,markershape=:star8)
plot!(avquatileHYsmSB,clquantileHYsmSB[1,:],fillrange=clquantileHYsmSB[3,:],fillalpha = 0.35, c=9, label="",lw = 0)
plot!(avquatileHYsmSB,clquantileHYsmSB[2,:],label=L"\mathrm{Hyperbolic \,\, }α=1.125, v=0.5, ζ=1 \mathrm{\,\, (SWSB)}",c=9,lw = 1,markershape=:star8)



    xlabel!(L"\mathrm{Average\,\, degree}")
    ylabel!(L"\mathrm{Temporal \,\, Clustering\,\, Coefficient}")
    title!(L"\mathrm{302\,\, nodes}")
    png("/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/optimizatedClustering300")


# #PATHLENGTH
    


p2=plot(avquatile,plquantile[1,:],fillrange=plquantile[3,:],fillalpha = 0.35, c = 1, label="",lw = 0,yaxis=:log,palette=:seaborn_colorblind)
plot!(avquatile,plquantile[2,:],label=L"\mathrm{Real\,\, data\,\, subjects\,\, Schaefer}",markershape=:circle,lw = 2,c=1)
plot!(avquatile,plquantileRE[1,:],fillrange=plquantileRE[3,:],fillalpha = 0.35, c=2, label="",lw = 0)
plot!(avquatile,plquantileRE[2,:],label=L"\mathrm{Random \,\, Edges}",markershape=:pentagon,lw = 2,c=2)
plot!(avquatile,plquantileRT[1,:],fillrange=plquantileRT[3,:],fillalpha = 0.35, c=3, label="",lw = 0)
plot!(avquatile,plquantileRT[2,:],label=L"\mathrm{Random \,\,Times}",markershape=:hexagon,lw = 2,c=3)

plot!(avquatileEU,plquantileEU[1,:],fillrange=plquantileEU[3,:],fillalpha = 0.35, c=4, label="",lw = 0)
plot!(avquatileEU,plquantileEU[2,:],label=L"\mathrm{Geometric \,\,on\,\, torus\,\, v=0.9}",markershape=:rect,lw = 2,c=4)
plot!(avquatileEUB,plquantileEUB[1,:],fillrange=plquantileEUB[3,:],fillalpha = 0.35, c=5, label="",lw = 0)
plot!(avquatileEUB,plquantileEUB[2,:],label=L"\mathrm{Geometric\,\,} v=0.9",markershape=:rect,lw = 2,c=5,legend=:topright,legendfontsize=5,xlims=(0, 165),grid=false,dpi=1200)
plot!(avquatileHYsm,plquantileHYsm[1,:],fillrange=plquantileHYsm[3,:],fillalpha = 0.35, c=10, label="",lw = 0)
plot!(avquatileHYsm,plquantileHYsm[2,:],label=L"\mathrm{Hyperbolic \,\, }α=1.120, v=0.6, ζ=1 \mathrm{\,\, (SW)}",c=10,lw = 1,markershape=:star8)
plot!(avquatileHYsmSB,plquantileHYsmSB[1,:],fillrange=plquantileHYsmSB[3,:],fillalpha = 0.35, c=9, label="",lw = 0)
plot!(avquatileHYsmSB,plquantileHYsmSB[2,:],label=L"\mathrm{Hyperbolic \,\, }α=1.125, v=0.5, ζ=1 \mathrm{\,\, (SWSB)}",c=9,lw = 1,markershape=:star8)


    xlabel!(L"\mathrm{Average\,\, degree}")
    ylabel!(L"\mathrm{Temporal\,\, Path \,\,Length}")
    title!(L"\mathrm{302\,\, nodes}")
    png("/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/optimizatedPathLenght300")


# #TEMPORALCORRELATIONCOEFFICIENT

p3=plot(avquatile,tccquantile[1,:],fillrange=tccquantile[3,:],fillalpha = 0.35, c = 1, label="",lw = 0,palette=:seaborn_colorblind)
plot!(avquatile,tccquantile[2,:],label=L"\mathrm{Real\,\, data\,\, subjects\,\, Schaefer}",markershape=:circle,lw = 2,c=1)
plot!(avquatile,tccquantileRE[1,:],fillrange=tccquantileRE[3,:],fillalpha = 0.35, c=2, label="",lw = 0)
plot!(avquatile,tccquantileRE[2,:],label=L"\mathrm{Random \,\, Edges}",markershape=:pentagon,lw = 2,c=2)
plot!(avquatile,tccquantileRT[1,:],fillrange=tccquantileRT[3,:],fillalpha = 0.35, c=3, label="",lw = 0)
plot!(avquatile,tccquantileRT[2,:],label=L"\mathrm{Random \,\,Times}",markershape=:hexagon,lw = 2,c=3)
plot!(avquatileEU,tccquantileEU[1,:],fillrange=tccquantileEU[3,:],fillalpha = 0.35, c=4, label="",lw = 0)
plot!(avquatileEU,tccquantileEU[2,:],label=L"\mathrm{Geometric \,\,on\,\, torus\,\, v=0.9}",markershape=:rect,lw = 2,c=4)
plot!(avquatileEUB,tccquantileEUB[1,:],fillrange=tccquantileEUB[3,:],fillalpha = 0.35, c=5, label="",lw = 0)
plot!(avquatileEUB,tccquantileEUB[2,:],label=L"\mathrm{Geometric\,\,} v=0.9",markershape=:rect,lw = 2,c=5,legend=:bottomright,legendfontsize=5,xlims=(0, 165),grid=false,dpi=1200)
plot!(avquatileHYsm,tccquantileHYsm[1,:],fillrange=tccquantileHYsm[3,:],fillalpha = 0.35, c=10, label="",lw = 0)
plot!(avquatileHYsm,tccquantileHYsm[2,:],label=L"\mathrm{Hyperbolic \,\, }α=1.120, v=0.6, ζ=1 \mathrm{\,\, (SW)}",c=10,lw = 1,markershape=:star8)
plot!(avquatileHYsmSB,tccquantileHYsmSB[1,:],fillrange=
tccquantileHYsmSB[3,:],fillalpha = 0.35, c=9, label="",lw = 0)
plot!(avquatileHYsmSB,tccquantileHYsmSB[2,:],label=L"\mathrm{Hyperbolic \,\, }α=1.125, v=0.5, ζ=1 \mathrm{\,\, (SWSB)}",c=9,lw = 1,markershape=:star8)


    xlabel!(L"\mathrm{Average\,\, degree}")
    ylabel!(L"\mathrm{Temporal \,\,Correlation \,\,Coefficient}")
    title!(L"\mathrm{302\,\, nodes}")
    png("/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/optimizatedTemporalCorrelationCoefficient300")

#SMALL WORLDNESS SIZEMORE AND BASSET
    
p4=plot(avquatile,smSBquantile[1,:],fillrange=smSBquantile[3,:],fillalpha = 0.35, c = 1, label="",lw = 0,palette=:seaborn_colorblind)
plot!(avquatile,smSBquantile[2,:],label=L"\mathrm{Real\,\, data\,\, subjects\,\, Schaefer}",markershape=:circle,lw = 2,c=1)
plot!(avquatile,smSBquantileRE[1,:],fillrange=smSBquantileRE[3,:],fillalpha = 0.35, c=2, label="",lw = 0)
plot!(avquatile,smSBquantileRE[2,:],label=L"\mathrm{Random \,\, Edges}",markershape=:pentagon,lw = 2,c=2)
plot!(avquatile,smSBquantileRT[1,:],fillrange=smSBquantileRT[3,:],fillalpha = 0.35, c=3, label="",lw = 0)
plot!(avquatile,smSBquantileRT[2,:],label=L"\mathrm{Random \,\,Times}",markershape=:hexagon,lw = 2,c=3)
plot!(avquatileEU,smSBquantileEU[1,:],fillrange=smSBquantileEU[3,:],fillalpha = 0.35, c=4, label="",lw = 0)
plot!(avquatileEU,smSBquantileEU[2,:],label=L"\mathrm{Geometric \,\,on\,\, torus\,\, v=0.9}",markershape=:rect,lw = 2,c=4)
plot!(avquatileEUB,smSBquantileEUB[1,:],fillrange=smSBquantileEUB[3,:],fillalpha = 0.35, c=5, label="",lw = 0)
plot!(avquatileEUB,smSBquantileEUB[2,:],label=L"\mathrm{Geometric\,\,} v=0.9",markershape=:rect,lw = 2,c=5,legend=:topleft,legendfontsize=5,xlims=(0, 165),grid=false,dpi=1200)
plot!(avquatileHYsm,smSBquantileHYsm[1,:],fillrange=smSBquantileHYsm[3,:],fillalpha = 0.35, c=10, label="",lw = 0)
plot!(avquatileHYsm,smSBquantileHYsm[2,:],label=L"\mathrm{Hyperbolic \,\, }α=1.120, v=0.6, ζ=1 \mathrm{\,\, (SW)}",c=10,lw = 1,markershape=:star8)
plot!(avquatileHYsmSB,smSBquantileHYsmSB[1,:],fillrange=smSBquantileHYsmSB[3,:],fillalpha = 0.35, c=9, label="",lw = 0)
plot!(avquatileHYsmSB,smSBquantileHYsmSB[2,:],label=L"\mathrm{Hyperbolic \,\, }α=1.125, v=0.5, ζ=1 \mathrm{\,\, (SWSB)}",c=9,lw = 1,markershape=:star8)
xlabel!(L"\mathrm{Average\,\, degree}")
ylabel!(L"\mathrm{Temporal \,\, Small \,\,Worldness\,\, S&B}")
title!(L"\mathrm{302 \,\,nodes}")
png("/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/optimizatedSmallWorldnessSizemoreBassett300")
