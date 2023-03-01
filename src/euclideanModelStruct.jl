using Random, Distributions,LinearAlgebra

abstract type TemporalNetwork end

struct EuclideanTemporalNetwork <: TemporalNetwork
    coordinates::Array{Float64, 3}
    adjacency::Array{Int64, 3}
    number_nodes::Int64
    number_snapshots::Int64
    ϵ::Float64
    d::Int64
    r::Float64

    function EuclideanTemporalNetwork(number_nodes::Int64,number_snapshots::Int64,ϵ::Float64,d::Int64)
        coordinates=Array{Float64}(undef, number_nodes,2,number_snapshots)
        adjacency=fill(0,(number_nodes,number_nodes,number_snapshots))
        coordinates=initializeNodes(coordinates,number_nodes)
        coordinates=updateNodes(coordinates,number_nodes,number_snapshots,ϵ)
        r=sqrt(d/(pi*(number_nodes-1))) 
        adjacency=fillAdjacencyMatrix(adjacency,coordinates,number_nodes,number_snapshots,r)
        new(coordinates,adjacency,number_nodes,number_snapshots,ϵ,d,r)
    end
end

function initializeNodes(coordinates::Array{Float64, 3},number_nodes::Int64)::Array{Float64, 3}
    for i in 1:number_nodes
        coordinates[i,:,1]=[rand(),rand()]
    end
    return coordinates 
end

function updateNodes(coordinates::Array{Float64, 3},number_nodes::Int64,number_snapshots::Int64,ϵ::Float64)::Array{Float64, 3}
    for t in 2:number_snapshots
        for i in 1:number_nodes
            r=rand(Uniform(-ϵ,ϵ))
            theta=rand(Uniform(0,2*pi))
            coordinates[i,1,t]=abs(coordinates[i,1,t-1]+r*cos(theta))%1
            coordinates[i,2,t]=abs(coordinates[i,2,t-1]+r*sin(theta))%1
        end
    end
    return coordinates
end

function distanceEU(nodeA::Array{Float64, 1},nodeB::Array{Float64, 1})::Float64
    dist=(1/2) .- abs.((1/2) .- abs.(nodeA-nodeB))
    return sqrt(dist[1]^2+dist[2]^2)
end

function fillAdjacencyMatrix(adjacency::Array{Int64, 3},coordinates::Array{Float64, 3},number_nodes::Int64,number_snapshots::Int64,r::Float64)::Array{Int64, 3}
    for t in 1:number_snapshots
        for j in 1:number_nodes
            for i in 1:number_nodes
                if distanceEU(coordinates[i,:,t],coordinates[j,:,t])<=r
                    adjacency[i,j,t]=1
                end
            end
        end
    end
    return adjacency
end




