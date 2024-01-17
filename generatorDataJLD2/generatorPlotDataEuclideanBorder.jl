using NPZ, Random, Distributions, LinearAlgebra, Permutations, Plots, JLD2, BenchmarkTools

include("../src/euclideanBorderModelStruct.jl")
include("../src/measuresStruct.jl")

radiusRange = collect(0:0.05:0.8)
velocities = collect(0.1:0.1:0.9)

Threads.@threads for n in 1:10
    avEU = zeros(length(velocities), length(radiusRange))
    clEU = zeros(length(velocities), length(radiusRange))
    plEU = zeros(length(velocities), length(radiusRange))
    tccEU = zeros(length(velocities), length(radiusRange))

    for i in 1:length(velocities)
        for j in 1:length(radiusRange)
            eu = EuclideanBorderTemporalNetwork(302, 27, velocities[i], radiusRange[j])
            avEU[i, j] = averageDegree(eu)
            plEU[i, j] = temporalPathLength(eu)
            clEU[i, j] = clusteringCoeff(eu)
            tccEU[i, j] = temporalCorrelationCoefficient(eu)
        end
    end

    @save "data/euclideanBorder_300_r_0_005_08_$(n).jld2" avEU clEU plEU tccEU radiusRange velocities
end