include("graphs.jl")
include("../src/measuresStruct.jl")

# plots to show the difference between the temporal clustering and the temporal correlation coefficient in two small examples

using Plots

ng = 30
number_snapshots = 50

graphs = Vector{StarGraph}(undef, ng)
tcc = Vector{Float64}(undef, ng)
cl = Vector{Float64}(undef, ng)

for i in 1:ng
    graphs[i] = StarGraph(rand(5:1:10), number_snapshots)
end


cl .= clusteringCoeff.(graphs)
tcc .= temporalCorrelationCoefficient.(graphs)


p1 = plot((1:ng), [cl, tcc], labels=["temporal clustering" "temporal correlation coefficient"], title="StarGraph with number_snapshots=$(number_snapshots)", xlabel="graphs", ylabel="value", legend=:bottomright, lw=2)


graphs = Vector{SwappingTriangle}(undef, ng)
tcc = Vector{Float64}(undef, ng)
cl = Vector{Float64}(undef, ng)

for i in 1:ng
    graphs[i] = SwappingTriangle(number_snapshots)
end


cl .= clusteringCoeff.(graphs)
tcc .= temporalCorrelationCoefficient.(graphs)


p2 = plot((1:ng), [cl, tcc], labels=["temporal clustering" "temporal correlation coefficient"], title="SwappingTriangle with number_snapshots=$(number_snapshots)", xlabel="graphs", ylabel="value", legend=:bottomright, lw=2)