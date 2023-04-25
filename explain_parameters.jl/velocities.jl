using NPZ, Random, Distributions,LinearAlgebra,Permutations,Plots,JLD2,BenchmarkTools

include("../src/hyperbolicModelStruct.jl")
include("../src/measuresStruct.jl")

#αrange=collect(0.5:0.025:1.2)
αrange=collect(0.5:0.025:1.4)
#Rrange=collect(0.5:0.5:14)
Rrange=collect(0:0.5:18)
velocities=collect(0.1:0.1:0.9)

res=zeros(length(velocities),26)
for k in 1:length(velocities)
    hy= HyperbolicTemporalNetwork(302,27,αrange[1],Rrange[10],velocities[k],1.0)
    for i in 1:26
        res[k,i]=sum(hy.adjacency[:,:,i].==hy.adjacency[:,:,i+1])/(302*302)*100
    end
end

points=mapslices(x->mean(x),res,dims=2)

p=plot(velocities,points,grid=false,xticks = ([ 0.2, 0.4, 0.6, 0.8], [ L"0.2", L"0.4", L"0.6", L"0.8"]),yticks = ([70, 75,80, 85, 90,95 ], [L"70\%", L"75\%", L"80\%", L"85\%", L"90\%",L"95\%"]),color=:black,lw=1.5,legend=false,xlabel=L"Velocity",ylabel=L"\mathrm{Percentage\,\, equal \,\, consecutive \,\,edges}", dpi=1200)
savefig(p, "/user/aurossi/home/mri_networks/TemporalBrainNetworksCode/images/velocities.png")
