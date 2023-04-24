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
            adjacency[:, :, t] .= Matrix(Graphs.adjacency_matrix(Graphs.watts_strogatz(number_nodes, degree, β)))
        end
        new(adjacency, number_nodes, number_snapshots, degree, β)
    end
end


function adj_watts_strogatz_temporal(
    n::Integer,
    k::Integer,
    β::Real;
    is_directed::Bool=false,
    remove_edges::Bool=true,
)
    @assert k < n
    adj = zeros(n, n, div(k, 2)*n+1)
    # If we have n - 1 neighbors (exactly k/2 on each side), then the graph is
    # necessarily complete. No need to run the Watts-Strogatz procedure:
    if k == n - 1 && iseven(k)
        return is_directed ? complete_digraph(n) : complete_graph(n)
    end

    g = is_directed ? SimpleDiGraph(n) : SimpleGraph(n)

    # The ith next vertex, in clockwise order.
    # (Reduce to zero-based indexing, so the modulo works, by subtracting 1
    # before and adding 1 after.)
    @inline target(s, i) = ((s + i - 1) % n) + 1

    # Phase 1: For each step size i, add an edge from each vertex s to the ith
    # next vertex, in clockwise order.
    t=2
    for i in 1:div(k, 2), s in 1:n
        add_edge!(g, s, target(s, i))
    end
    adj[:,:,1] .= Matrix(Graphs.adjacency_matrix(g))
    for i in 1:div(k, 2), s in 1:n

        # We only rewire with a probability β, and we only worry about rewiring
        # if there is some vertex not connected to s; otherwise, the only valid
        # rewiring is to reconnect to the ith next vertex, and there is no work
        # to do.
        (rand() < β && degree(g, s) < n - 1) || continue

        t = target(s, i)

        while true
            d = rand(1:n)          # Tentative new target
            d == s && continue          # Self-loops prohibited
            d == t && break             # Rewired to original target
            if add_edge!(g, s, d)       # Was this valid (i.e., unconnected)?
                remove_edges && rem_edge!(g, s, t)     # True rewiring: Delete original edge
                break                                   # We found a valid target
            end
        end
        adj[:,:,t] .= Matrix(Graphs.adjacency_matrix(g))
        t+=1
    end
    return adj[:,:,1:t-1]
end


struct WattsStrogatzTemporal <: TemporalNetwork
    adjacency::Array{Int64,3}
    number_nodes::Int64
    number_snapshots::Int64
    degree::Int64
    β::Float64


    function WattsStrogatzTemporal(number_nodes::Int64, number_snapshots::Int64, degree::Int64, β::Float64)
        adj=adj_watts_strogatz_temporal(number_nodes, degree, β)
        new(adj, number_nodes,size(adj,3) , degree, β)
    end
end
