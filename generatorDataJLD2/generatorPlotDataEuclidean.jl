using NPZ, Random, Distributions, LinearAlgebra, Permutations, Plots, JLD2, BenchmarkTools, DelimitedFiles

include("../src/euclideanModelStruct.jl")
include("../src/measuresStruct.jl")

thresholds = append!(collect(0.2:0.05:0.90), collect(0.92:0.02:0.98))
velocities = collect(0.1:0.1:0.9)
subjects = readdlm("../src/filtered-subjects-mod.txt", Int)

avall = zeros(length(subjects), size(thresholds)[1])

for i in 1:length(subjects)
    @load "/data/Hyperbrain/$(subjects[i])/$(subjects[i])_schaefer_300_thre_02_005_098_OR_RT_RE_ws60_wo30.jld2" av cl pl tcc

    avall[i, :] = av
end

Threads.@threads for n in 1:10
    avEU = zeros(length(velocities), length(thresholds))
    clEU = zeros(length(velocities), length(thresholds))
    plEU = zeros(length(velocities), length(thresholds))
    tccEU = zeros(length(velocities), length(thresholds))

    for i in 1:length(velocities)
        for j in 1:length(thresholds)
            eu = EuclideanTemporalNetwork(302, 27, velocities[i], round(Int, mean(avall[:, j])))
            avEU[i, j] = averageDegree(eu)
            plEU[i, j] = temporalPathLength(eu)
            clEU[i, j] = clusteringCoeff(eu)
            tccEU[i, j] = temporalCorrelationCoefficient(eu)
        end
    end
    @save "../data/euclidean_300_thre_02_005_098_$(n).jld2" avEU clEU plEU tccEU thresholds velocities

end