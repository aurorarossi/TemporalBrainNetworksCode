using NPZ, Random, Distributions,LinearAlgebra,Permutations,Plots,JLD2,BenchmarkTools

include("../src/hyperbolicModelStruct.jl")
include("../src/measuresStruct.jl")

#αrange=collect(0.5:0.025:1.2)
αrange=collect(0.5:0.025:1.4)
#Rrange=collect(0.5:0.5:14)
Rrange=collect(0:0.5:18)
velocities=collect(0.1:0.1:0.9)
#Threads.@threads for n in 1:10
    tccHY=zeros(length(velocities),length(Rrange),length(αrange))
    avHY=zeros(length(velocities),length(Rrange),length(αrange))
    clHY=zeros(length(velocities),length(Rrange),length(αrange))
    plHY=zeros(length(velocities),length(Rrange),length(αrange))
    for a in 1:length(αrange)
        @time begin
            for i in 1:length(velocities)
                for j in 1:length(Rrange)
                    hy= HyperbolicTemporalNetwork(302,27,αrange[a],Rrange[j],velocities[i],1.0)
                    tccHY[i,j,a]=temporalCorrelationCoefficient(hy)
                    avHY[i,j,a]=averageDegree(hy)
                    plHY[i,j,a]=temporalPathLength(hy)
                    clHY[i,j,a]=clusteringCoeff(hy)
                end
                println("v:$(i)")
            end

            #@save "hyperbolic_300_alpha_$(round(Int,αrange[a]*100))_R_05_05_7_v_01_01_09_$(n).jld2" avHY clHY plHY tccHY αrange Rrange velocities
            println("alpha:$(a)")
        end
    end
@save "/user/aurossi/home/mri_networks/TemporalBrainNetworksCode/data/hyperbolic_300.jld2" avHY clHY plHY tccHY Rrange αrange velocities