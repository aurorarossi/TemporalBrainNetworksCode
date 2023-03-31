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
degree =  75

betas = [0.0001,0.0005,0.001,0.01,0.02,0.07,0.1,0.2,0.5,1.0]

graphs =  Vector{WattsStrogatz}(undef, 10)
tcc = Vector{Float64}(undef, 10)
cl = Vector{Float64}(undef, 10)
cl_graphs = Vector{Float64}(undef, 10)
for i in 1:10
    graphs[i] = WattsStrogatz(number_nodes, number_snapshots, degree , i/10)
end

for i in 1:10
    cl[i] = clusteringCoeff(graphs[i])
    tcc[i] = temporalCorrelationCoefficient(graphs[i])
    #cl_graphs[i] = clusteringCoeffGraphs(graphs[i])/clusteringCoeffGraphs(graph_β)
end

plot([cl, tcc],labels=["temporal clustering " "temporal correlation coefficient"],title="number_snapshots=$(number_snapshots)",xlabel="β",ylabel="value",legend=:topright)