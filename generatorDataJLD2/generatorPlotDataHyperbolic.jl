using NPZ, Random, Distributions,LinearAlgebra,Permutations,Plots,JLD2,BenchmarkTools

include("../src/hyperbolicModelStruct.jl")
include("../src/measuresStruct.jl")

#αrange=collect(0.5:0.025:1.2)
αrange=collect(0.5:0.025:1.4)
#Rrange=collect(0.5:0.5:14)
Rrange=collect(0:0.5:18)
velocities=collect(0.1:0.1:0.9)
velocities=[0.8]
αrange=[1.05]
Threads.@threads for n in 1:10
    for a in 1:length(αrange)
        @time begin
            tccHY=zeros(length(velocities),length(Rrange))
            avHY=zeros(length(velocities),length(Rrange))
            clHY=zeros(length(velocities),length(Rrange))
            plHY=zeros(length(velocities),length(Rrange))
            for i in 1:length(velocities)
                for j in 1:length(Rrange)
                    hy= HyperbolicTemporalNetwork(302,27,αrange[a],Rrange[j],velocities[i],1.0)
                    tccHY[i,j]=temporalCorrelationCoefficient(hy)
                    avHY[i,j]=averageDegree(hy)
                    plHY[i,j]=temporalPathLength(hy)
                    clHY[i,j]=clusteringCoeff(hy)
                end
                println("v:$(i)")
            end
            @save "/user/aurossi/home/mri_networks/TemporalBrainNetworksCode/data/hyperbolic_300_alpha_$(round(Int,αrange[a]*100))_R_0_05_18_v_$(round(Int,velocities[1]*100))_$(n).jld2" avHY clHY plHY tccHY Rrange αrange velocities
            println("alpha:$(a)")
        end
    end
end