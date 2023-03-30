using Graphs

abstract type TemporalNetwork end

struct WattsStrogatz <: TemporalNetwork
    adjacency::Array{Int64,3}
    number_nodes::Int64
    number_snapshots::Int64
    degree::Int64
    β::Float64


    function WattsStrogatz(number_nodes::Int64, number_snapshots::Int64, degree::Int64, β::Float64)
        adjacency=zeros(number_nodes,number_nodes,number_snapshots)
        for t in 1:number_snapshots
            adjacency[:, :, t] .= Matrix(adjacency_matrix(watts_strogatz(number_nodes, degree, β)))
        end
        new(adjacency, number_nodes, number_snapshots, degree, β)
    end
end