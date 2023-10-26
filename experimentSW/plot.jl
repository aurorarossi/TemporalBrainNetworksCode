include("star_graph.jl")
include("../src/measuresStruct.jl")
using Plots
ng=30
number_snapshots =50
graphs =  Vector{StarGraph}(undef, ng)
tcc = Vector{Float64}(undef, ng)
cl = Vector{Float64}(undef, ng)
cl_graphs = Vector{Float64}(undef, ng)
for i in 1:ng
    graphs[i] = StarGraph(rand(5:1:10),number_snapshots )
end

for i in 1:ng
    cl[i] = clusteringCoeff(graphs[i])
    tcc[i] = temporalCorrelationCoefficient(graphs[i])
    #cl_graphs[i] = clusteringCoeffGraphs(graphs[i])/clusteringCoeffGraphs(graph_β)
end

p1=plot((1:ng),[cl, tcc],labels=["temporal clustering " "temporal correlation coefficient"],title="number_snapshots=$(number_snapshots)",xlabel="β",ylabel="value",legend=:bottomright)


graphs =  Vector{SwappingTriangle}(undef, ng)
tcc = Vector{Float64}(undef, ng)
cl = Vector{Float64}(undef, ng)
cl_graphs = Vector{Float64}(undef, ng)
for i in 1:ng
    graphs[i] = SwappingTriangle(number_snapshots)
end

for i in 1:ng
    cl[i] = clusteringCoeff(graphs[i])
    tcc[i] = temporalCorrelationCoefficient(graphs[i])
    #cl_graphs[i] = clusteringCoeffGraphs(graphs[i])/clusteringCoeffGraphs(graph_β)
end

p2=plot((1:ng),[cl, tcc],labels=["temporal clustering " "temporal correlation coefficient"],title="number_snapshots=$(number_snapshots)",xlabel="β",ylabel="value",legend=:bottomright)