using NPZ, Random, Distributions, LinearAlgebra, Permutations, Plots, JLD2, DelimitedFiles

include("../src/permutatedModelsStruct.jl")
include("../src/measuresStruct.jl")

abstract type TemporalNetwork end

struct OriginalTemporalNetwork <: TemporalNetwork
    adjacency::Array{Int64,3}
    number_nodes::Int64
    number_snapshots::Int64
end

thresholds = append!(collect(0.2:0.05:0.90), collect(0.92:0.02:0.98))
subjects = readdlm("../src/filtered-subjects.txt", Int)

for subject in subjects
    av = zeros(size(thresholds)[1])
    avRT = zeros(size(thresholds)[1])
    avRE = zeros(size(thresholds)[1])

    cl = zeros(size(thresholds)[1])
    clRT = zeros(size(thresholds)[1])
    clRE = zeros(size(thresholds)[1])

    pl = zeros(size(thresholds)[1])
    plRT = zeros(size(thresholds)[1])
    plRE = zeros(size(thresholds)[1])

    tcc = zeros(size(thresholds)[1])
    tccRT = zeros(size(thresholds)[1])
    tccRE = zeros(size(thresholds)[1])


    Threads.@threads for i in 1:size(thresholds)[1]
        path = "/data/Hyperbrain/$(subject)/$(subject)_ws60_wo30_t0_pearson_0_schaefer_300_b$(trunc(Int, (thresholds[i]*100))).npy"
        if isfile(path)
            network = npzread(path)
            networkpermuted = round.(Int, permutedims(network, (2, 3, 1)))
            networkStruct = OriginalTemporalNetwork(networkpermuted, size(networkpermuted)[1], size(networkpermuted)[3])

            av[i] = averageDegree(networkStruct)
            cl[i] = clusteringCoeff(networkStruct)
            pl[i] = temporalPathLength(networkStruct)
            tcc[i] = temporalCorrelationCoefficient(networkStruct)

            RT = RTTemporalNetwork(networkpermuted)
            avRT[i] = averageDegree(RT)
            clRT[i] = clusteringCoeff(RT)
            plRT[i] = temporalPathLength(RT)
            tccRT[i] = temporalCorrelationCoefficient(RT)

            RE = RETemporalNetwork(networkpermuted)
            avRE[i] = averageDegree(RE)
            clRE[i] = clusteringCoeff(RE)
            plRE[i] = temporalPathLength(RE)
            tccRE[i] = temporalCorrelationCoefficient(RE)

        end
    end
    @save "/data/$(subject)/$(subject)_schaefer_300_thre_02_005_098_OR_RT_RE_ws60_wo30.jld2" av cl pl tcc avRT clRT plRT tccRT avRE clRE plRE tccRE
end


