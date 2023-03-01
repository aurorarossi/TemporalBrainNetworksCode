using NPZ, Random, Distributions,LinearAlgebra,Permutations,Plots,JLD2,Statistics, DelimitedFiles,Dierckx

thresholds=append!(collect(0.2:0.05:0.90),collect(0.92:0.02:0.98))
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

p=plot()

for i in 1:length(subjects)
    @load "/data/Hyperbrain/$(subjects[i])/$(subjects[i])_schaefer_300_thre_02_005_098_OR_RT_RE_ws60_wo30.jld2" av cl pl tcc avRT clRT plRT tccRT avRE clRE plRE tccRE
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

varvect=clquantileEU[2,:]./plquantileEU[2,:]
spline=Spline1D(sort(reverse(avquatileEU[:])),reverse(varvect[:]))
splineorg=Spline1D(sort(reverse(avquatile[:])),reverse(smquantile[2,:]))
range=collect(0:1:165)
evspline=evaluate(spline,range)
evsplineorg=evaluate(splineorg,range)
diffspline=Spline1D(range,abs.(evspline-evsplineorg)) #HEREE
alld=integrate(diffspline,0,165)
println("Distance smallworldness EU TORUS")
println(alld)
plot!(range,evspline)
plot!(range,evsplineorg)

varvect=clquantileEUB[2,:]./plquantileEUB[2,:]
spline=Spline1D(sort((avquatileEUB[:])),(varvect[:]))
splineorg=Spline1D(sort(reverse(avquatile[:])),reverse(smquantile[2,:]))
range=collect(0:1:165)
evspline=evaluate(spline,range)
evsplineorg=evaluate(splineorg,range)
diffspline=Spline1D(range,abs.(evspline-evsplineorg)) #HEREE
alld=integrate(diffspline,0,165)
println("Distance smallworldness EU SQUARE")
println(alld)
plot!(range,evspline)
plot!(range,evsplineorg)

varvect=clquantileRT[2,:]./plquantileRT[2,:]   
spline=Spline1D(sort(reverse(avquatileRT[:])),reverse(varvect[:]))
splineorg=Spline1D(sort(reverse(avquatile[:])),reverse(smquantile[2,:]))
range=collect(0:1:165)
evspline=evaluate(spline,range)
evsplineorg=evaluate(splineorg,range)
diffspline=Spline1D(range,abs.(evspline-evsplineorg)) #HEREE
alld=integrate(diffspline,0,165)
println("Distance smallworldness RT")
println(alld)
plot!(range,evspline)
plot!(range,evsplineorg)

varvect=clquantileRE[2,:]./plquantileRE[2,:]   
spline=Spline1D(sort(reverse(avquatileRE[:])),reverse(varvect[:]))
splineorg=Spline1D(sort(reverse(avquatile[:])),reverse(smquantile[2,:]))
range=collect(0:1:165)
evspline=evaluate(spline,range)
evsplineorg=evaluate(splineorg,range)
diffspline=Spline1D(range,abs.(evspline-evsplineorg)) #HEREE
alld=integrate(diffspline,0,165)
println("Distance smallworldness RE")
println(alld)
plot!(range,evspline)
plot!(range,evsplineorg)
# αrange=collect(0.5:0.025:1.2)
 velocities=collect(0.1:0.1:0.9)

αrange=[0.65]
a=1
indexv=6

@load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/hyperbolic_300_alpha_$(round(Int,αrange[a]*100))_R_05_05_14_v_01_01_09_correct.jld2" avHY clHY plHY tccHY 
avHYct=zeros(length(velocities),length(thresholds))
varvect=zeros(length(velocities),length(thresholds))
avHYvect=avHY[:,:]
varvect=clHY[:,:]./plHY[:,:]   #HERE
avHYvect=avHYvect[indexv,:]
varvect=varvect[indexv,:]

spline=Spline1D(sort(reverse(avHYvect[:])),reverse(varvect[:]))
splineorg=Spline1D(sort(reverse(avquatile[:])),reverse(smquantile[2,:]))
range=collect(0:1:165)
evspline=evaluate(spline,range)
evsplineorg=evaluate(splineorg,range)
diffspline=Spline1D(range,abs.(evspline-evsplineorg)) #HEREE
alld=integrate(diffspline,0,165)
plot!(range,evspline)
plot!(range,evsplineorg)



println("Distance smallworldness hype sm")
println(alld)

αrange=[0.975]
a=1
indexv=9

@load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/hyperbolic_300_alpha_$(round(Int,αrange[a]*100))_R_05_05_14_v_01_01_09_correct.jld2" avHY clHY plHY tccHY 
avHYct=zeros(length(velocities),length(thresholds))
varvect=zeros(length(velocities),length(thresholds))
avHYvect=avHY[:,:]
varvect=clHY[:,:]./plHY[:,:]   #HERE
avHYvect=avHYvect[indexv,:]
varvect=varvect[indexv,:]

spline=Spline1D(sort(reverse(avHYvect[:])),reverse(varvect[:]))
splineorg=Spline1D(sort(reverse(avquatile[:])),reverse(smquantile[2,:]))
range=collect(0:1:165)
evspline=evaluate(spline,range)
evsplineorg=evaluate(splineorg,range)
diffspline=Spline1D(range,abs.(evspline-evsplineorg)) #HEREE
alld=integrate(diffspline,0,165)
plot!(range,evspline)
plot!(range,evsplineorg)



println("Distance smallworldness hype SB")
println(alld)

################################################################
println()
println()
println()
println()
println()
println()

varvect=tccquantileEU[2,:]./plquantileEU[2,:]
spline=Spline1D(sort(reverse(avquatileEU[:])),reverse(varvect[:]))
splineorg=Spline1D(sort(reverse(avquatile[:])),reverse(smSBquantile[2,:]))
range=collect(0:1:165)
evspline=evaluate(spline,range)
evsplineorg=evaluate(splineorg,range)
diffspline=Spline1D(range,abs.(evspline-evsplineorg)) #HEREE
alld=integrate(diffspline,0,165)
println("Distance smallworldness Sizemore and Bassett EU TORUS")
println(alld)
plot!(range,evspline)
plot!(range,evsplineorg)

varvect=tccquantileEUB[2,:]./plquantileEUB[2,:]
spline=Spline1D(sort((avquatileEUB[:])),(varvect[:]))
splineorg=Spline1D(sort(reverse(avquatile[:])),reverse(smSBquantile[2,:]))
range=collect(0:1:165)
evspline=evaluate(spline,range)
evsplineorg=evaluate(splineorg,range)
diffspline=Spline1D(range,abs.(evspline-evsplineorg)) #HEREE
alld=integrate(diffspline,0,165)
println("Distance smallworldness Sizemore and Bassett EU SQUARE")
println(alld)
plot!(range,evspline)
plot!(range,evsplineorg)


varvect=tccquantileRT[2,:]./plquantileRT[2,:]   
spline=Spline1D(sort(reverse(avquatileRT[:])),reverse(varvect[:]))
splineorg=Spline1D(sort(reverse(avquatile[:])),reverse(smSBquantile[2,:]))
range=collect(0:1:165)
evspline=evaluate(spline,range)
evsplineorg=evaluate(splineorg,range)
diffspline=Spline1D(range,abs.(evspline-evsplineorg)) #HEREE
alld=integrate(diffspline,0,165)
println("Distance Sizemore and Bassett smallworldness RT")
println(alld)


varvect=tccquantileRE[2,:]./plquantileRE[2,:]   
spline=Spline1D(sort(reverse(avquatileRE[:])),reverse(varvect[:]))
splineorg=Spline1D(sort(reverse(avquatile[:])),reverse(smSBquantile[2,:]))
range=collect(0:1:165)
evspline=evaluate(spline,range)
evsplineorg=evaluate(splineorg,range)
diffspline=Spline1D(range,abs.(evspline-evsplineorg)) #HEREE
alld=integrate(diffspline,0,165)
println("Distance Sizemore and Bassett smallworldness RE")
println(alld)



αrange=[0.65]
a=1
indexv=6

@load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/hyperbolic_300_alpha_$(round(Int,αrange[a]*100))_R_05_05_14_v_01_01_09_correct.jld2" avHY clHY plHY tccHY 
avHYct=zeros(length(velocities),length(thresholds))
varvect=zeros(length(velocities),length(thresholds))
avHYvect=avHY[:,:]
varvect=tccHY[:,:]./plHY[:,:]   #HERE
avHYvect=avHYvect[indexv,:]
varvect=varvect[indexv,:]

spline=Spline1D(sort(reverse(avHYvect[:])),reverse(varvect[:]))
splineorg=Spline1D(sort(reverse(avquatile[:])),reverse(smSBquantile[2,:]))
range=collect(0:1:165)
evspline=evaluate(spline,range)
evsplineorg=evaluate(splineorg,range)
diffspline=Spline1D(range,abs.(evspline-evsplineorg)) #HEREE
alld=integrate(diffspline,0,165)
 



println("Distance Sizemore and Bassett smallworldness hype sm")
println(alld)


αrange=[0.975]
a=1
indexv=9

@load "/user/aurossi/home/mri_networks/hyperbolicbrains/optimization/hyperbolic_300_alpha_$(round(Int,αrange[a]*100))_R_05_05_14_v_01_01_09_correct.jld2" avHY clHY plHY tccHY 
avHYct=zeros(length(velocities),length(thresholds))
varvect=zeros(length(velocities),length(thresholds))
avHYvect=avHY[:,:]
varvect=tccHY[:,:]./plHY[:,:]   #HERE
avHYvect=avHYvect[indexv,:]
varvect=varvect[indexv,:]

spline=Spline1D(sort(reverse(avHYvect[:])),reverse(varvect[:]))
splineorg=Spline1D(sort(reverse(avquatile[:])),reverse(smSBquantile[2,:]))
range=collect(0:1:165)
evspline=evaluate(spline,range)
evsplineorg=evaluate(splineorg,range)
diffspline=Spline1D(range,abs.(evspline-evsplineorg)) #HEREE
alld=integrate(diffspline,0,165)
 



println("Distance Sizemore and Bassett smallworldness hype SB")
println(alld)







































