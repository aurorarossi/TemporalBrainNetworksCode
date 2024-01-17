abstract type TemporalNetwork end

function averageDegree(network::TemporalNetwork)::Float64
    averagedegree = 0
    for t in 1:network.number_snapshots
        averagesnap = 0
        for j in 1:network.number_nodes
            for i in 1:network.number_nodes
                averagesnap += network.adjacency[i, j, t]
            end
        end
        averagedegree += (averagesnap) / network.number_nodes - 1
    end
    return averagedegree / network.number_snapshots
end

function clustringCoeffOneSnap(t::Int64, network::TemporalNetwork)::Float64
    s = 0
    for i in 1:network.number_nodes
        for j in 1:network.number_nodes
            if i != j
                for k in 1:network.number_nodes
                    if j != k && k != i
                        s += network.adjacency[i, j, t] * network.adjacency[j, k, t] * network.adjacency[k, i, t]
                    end
                end
            end
        end
    end
    l = 0
    ki = 0
    for j in 1:network.number_nodes
        for i in 1:network.number_nodes
            i != j && (ki += network.adjacency[i, j, t])
        end
        l += ki * (ki - 1)
        ki = 0
    end
    gc = s / l
    l == 0 && (gc = 0.0)
    return gc
end

function clusteringCoeff(network::TemporalNetwork)::Float64
    s = 0
    for t in 1:network.number_snapshots
        s += clustringCoeffOneSnap(t, network)
    end
    return s / network.number_snapshots
end

function temporalCorrelationCoefficient(network::TemporalNetwork)::Float64
    ci = 0
    sum_ci = 0
    for i in 1:network.number_nodes
        st = 0
        for t in 1:network.number_snapshots-1
            n = 0
            for j in 1:network.number_nodes
                n += (network.adjacency[i, j, t] * network.adjacency[i, j, t+1])
            end
            d = sqrt(sum(network.adjacency[i, :, t]) * sum(network.adjacency[i, :, t+1]))
            if d != 0
                st += n / d
            end
        end
        ci = st / (network.number_snapshots - 1)
        sum_ci += ci
    end
    return sum_ci / network.number_nodes
end

function BFS(i::Int64, network::TemporalNetwork)::Array{Float64,1}
    visited = zeros(network.number_nodes)
    visited[i] = 1
    distance = network.number_snapshots .* ones(network.number_nodes)
    distance[i] = 0
    for j in 1:network.number_nodes
        if visited[j] == 0 && network.adjacency[i, j, 1] == 1
            distance[j] = 1
            visited[j] = 1
        end
    end
    for t in 2:network.number_snapshots
        for i in 1:length(visited)
            if visited[i] == 1
                for j in 1:network.number_nodes
                    if visited[j] == 0 && network.adjacency[i, j, t] == 1
                        if t < distance[j]
                            distance[j] = t
                            visited[j] = 1
                        end
                    end
                end
            end
        end
    end
    return distance
end

function meanL(distances::Array{Float64,2})::Float64
    s = 0
    for j in 1:size(distances)[2]
        for i in 1:size(distances)[1]
            if i != j
                s += distances[i, j]
            end
        end
    end
    return (1 / (size(distances)[1] * (size(distances)[1] - 1))) * s
end

function temporalPathLength(network::TemporalNetwork)::Float64
    matrixBFS = zeros(network.number_nodes, network.number_nodes)
    Threads.@threads for i in 1:network.number_nodes
        matrixBFS[:, i] = BFS(i, network)
    end
    meanL(matrixBFS)
end