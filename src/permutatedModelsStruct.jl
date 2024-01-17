using Random
abstract type TemporalNetwork end

struct RTTemporalNetwork <: TemporalNetwork
    adjacency::Array{Int64,3}
    number_nodes::Int64
    number_snapshots::Int64


    function RTTemporalNetwork(original_adjacency::Array{Int64,3})
        number_nodes = size(original_adjacency)[1]
        number_snapshots = size(original_adjacency)[3]
        adjacencyRT = zeros(number_nodes, number_nodes, number_snapshots)
        for i in 1:number_nodes
            for j in 1:i
                tpermuted = array(RandomPermutation(number_snapshots))
                for t in 1:number_snapshots
                    adjacencyRT[i, j, t] = adjacencyRT[j, i, t] = original_adjacency[i, j, tpermuted[t]]
                end
            end
        end
        new(adjacencyRT, number_nodes, number_snapshots)
    end
end

struct RETemporalNetwork <: TemporalNetwork
    adjacency::Array{Int64,3}
    number_nodes::Int64
    number_snapshots::Int64


    function RETemporalNetwork(original_adjacency::Array{Int64,3})
        number_nodes = size(original_adjacency)[1]
        number_snapshots = size(original_adjacency)[3]
        adjacencyRE = zeros(number_nodes, number_nodes, number_snapshots)
        for t in 1:number_snapshots
            for i in 1:number_nodes
                adjacencyRE[i, i, t] = 1
                for j in i+1:number_nodes
                    if original_adjacency[i, j, t] == 1
                        indexrandom = rand(1:number_nodes)
                        if rand() > 0.5
                            if (adjacencyRE[i, indexrandom, t] != 1) && indexrandom != i
                                adjacencyRE[i, indexrandom, t] = adjacencyRE[indexrandom, i, t] = 1
                            else
                                adjacencyRE[i, j, t] = adjacencyRE[j, i, t] = 1
                            end
                        else
                            if (adjacencyRE[indexrandom, j, t] != 1) && indexrandom != j
                                adjacencyRE[indexrandom, j, t] = adjacencyRE[j, indexrandom, t] = 1
                            else
                                adjacencyRE[i, j, t] = adjacencyRE[j, i, t] = 1
                            end

                        end
                    end
                end
            end
        end
        new(adjacencyRE, number_nodes, number_snapshots)
    end
end