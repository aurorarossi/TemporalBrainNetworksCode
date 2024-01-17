using Graphs

# struct for StarGraph, RandomDense and SwappingTriangle graphs

abstract type TemporalNetwork end

struct StarGraph <: TemporalNetwork
    adjacency::Array{Int64,3}
    number_nodes::Int64
    number_snapshots::Int64


    function StarGraph(number_nodes::Int64, number_snapshots::Int64)
        adjacency=zeros(number_nodes,number_nodes,number_snapshots)
        for t in 1:number_snapshots
            adjacency[:, :, t] .= Matrix(Graphs.adjacency_matrix(Graphs.star_graph(number_nodes)))
        end
        new(adjacency, number_nodes, number_snapshots)
    end
end

struct RandomDense <: TemporalNetwork
    adjacency::Array{Int64,3}
    number_nodes::Int64
    number_snapshots::Int64
    p::Float64


    function RandomDense(number_nodes::Int64, number_snapshots::Int64,p)
        adjacency=zeros(number_nodes,number_nodes,number_snapshots)
        for t in 1:number_snapshots
            adjacency[:, :, t] .= Matrix(Graphs.adjacency_matrix(Graphs.erdos_renyi(number_nodes,p)))
        end
        new(adjacency, number_nodes, number_snapshots)
    end
end

struct SwappingTriangle <: TemporalNetwork
    adjacency::Array{Int64,3}
    number_snapshots::Int64
    number_nodes::Int64


    function SwappingTriangle(number_snapshots::Int64)
        adjacency=zeros(6,6,number_snapshots)
        for t in 1:number_snapshots
            if t % 2 == 0
                adjacency[:, :, t] .= [0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0;0 0 0 0 0 0]
            else
                adjacency[:, :, t] .= [0 0 0 0 0 0; 0 0 0 0 0 0;0 0 0 0 0 0; 0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 0 0 0]
            end
        end
        new(adjacency, number_snapshots, 6)
    end
end