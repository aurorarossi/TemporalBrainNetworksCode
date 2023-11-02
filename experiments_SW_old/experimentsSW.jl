include("WattsStrogatz.jl")
include("../src/measuresStruct.jl")
using Plots

function clusteringCoeffGraphs(network::TemporalNetwork)::Float64
    s=0
    for t in 1:network.number_snapshots
        g = Graphs.Graph(network.adjacency[:,:,t])
        s+=Graphs.global_clustering_coefficient(g)
    end
    return s/network.number_snapshots
end

number_nodes = 150
number_snapshots = 50
deg =  20



graphs =  Vector{WattsStrogatz}(undef, 10)
tcc = Vector{Float64}(undef, 10)
cl = Vector{Float64}(undef, 10)
cl_graphs = Vector{Float64}(undef, 10)
for i in 1:10
    graphs[i] = WattsStrogatz(number_nodes, number_snapshots, deg , i/10)
end

for i in 1:10
    cl[i] = clusteringCoeff(graphs[i])
    tcc[i] = temporalCorrelationCoefficient(graphs[i])
    #cl_graphs[i] = clusteringCoeffGraphs(graphs[i])/clusteringCoeffGraphs(graph_β)
end

p1=plot((1:10)./10,[cl, tcc],labels=["temporal clustering " "temporal correlation coefficient"],title="number_snapshots=$(number_snapshots)",xlabel="β",ylabel="value",legend=:topright)

#Temporal

graphs =  Vector{WattsStrogatzTemporal}(undef, 10)
tcc = Vector{Float64}(undef, 10)
cl = Vector{Float64}(undef, 10)
cl_graphs = Vector{Float64}(undef, 10)
for i in 1:10
    graphs[i] = WattsStrogatzTemporal(number_nodes, number_snapshots, deg , i/10)
end

for i in 1:10
    cl[i] = clusteringCoeff(graphs[i])
    tcc[i] = temporalCorrelationCoefficient(graphs[i])
    #cl_graphs[i] = clusteringCoeffGraphs(graphs[i])/clusteringCoeffGraphs(graph_β)
end

p2=plot((1:10)./10,[cl, tcc],labels=["temporal clustering " "temporal correlation coefficient"],title="number_snapshots=$(number_snapshots)",xlabel="β",ylabel="value",legend=:bottomright)