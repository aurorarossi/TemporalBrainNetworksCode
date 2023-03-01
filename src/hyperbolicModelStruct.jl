using Random, Distributions

abstract type TemporalNetwork end

struct HyperbolicTemporalNetwork <: TemporalNetwork
    coordinates::Array{Float64, 3}
    probabilities::Array{Float64, 2}
    adjacency::Array{Int64, 3}
    number_nodes::Int64
    number_snapshots::Int64
    α::Float64
    R::Float64
    ϵ::Float64
    ζ::Float64

    function HyperbolicTemporalNetwork(number_nodes::Int64,number_snapshots::Int64,α::Float64,R::Float64,ϵ::Float64,ζ::Float64)
        coordinates=Array{Float64}(undef, number_nodes,2,number_snapshots)
        probabilities=Array{Float64}(undef, number_nodes,number_snapshots)
        adjacency=fill(0,(number_nodes,number_nodes,number_snapshots))

        coordinates,probabilities=initializeNodes(coordinates,probabilities,number_nodes,α,R)
        coordinates,probabilities=updateNodes(coordinates,probabilities,number_nodes,number_snapshots,α,R,ϵ)
        adjacency=fillAdjacencyMatrix(adjacency,coordinates,number_nodes,number_snapshots,R,ζ)
        new(coordinates,probabilities,adjacency,number_nodes,number_snapshots,α,R,ϵ,ζ)
    end
end

"""
    distanceHY(nodeA::Array{Float64, 1},nodeB::Array{Float64, 1}, ζ=1)::Float64

Compute the hyperbolic distance between two points.
"""
function distanceHY(nodeA::Array{Float64, 1},nodeB::Array{Float64, 1}, ζ=1)::Float64
    if nodeA!=nodeB
        d=acosh(cosh(ζ*nodeA[1])*cosh(ζ*nodeB[1])-(sinh(ζ*nodeA[1])*sinh(ζ*nodeB[1])*cos(pi-abs(pi-abs(nodeA[2]-nodeB[2])))))/ζ
        #println(cosh(BigFloat(ζ*nodeA[1]))==sinh(BigFloat(ζ*nodeA[1])))
    else
        d=0.0
    end
    return d
end

"""
    initializeNodes(coordinates::Array{Float64, 3},probabilities::Array{Float64,3},number_nodes::Int64,α::Float64,R::Float64)::Tuple{Array{Float64, 3},Array{Float64, 3}}

Initialize the points on the hyperboloid.
"""
function initializeNodes(coordinates::Array{Float64, 3},probabilities::Array{Float64,2},number_nodes::Int64,α::Float64,R::Float64)::Tuple{Array{Float64, 3},Array{Float64, 2}}
    for i in 1:number_nodes
        probabilities[i,1]=rand(Uniform(0,1))
        coordinates[i,1,1]=(1/α)*acosh(1+(cosh(α*R)-1)*probabilities[i,1])
        coordinates[i,2,1]=rand(Uniform(0,2*pi))
    end
    return coordinates, probabilities
end

"""
    updateNodes(coordinates::Array{Float64, 3},probabilities::Array{Float64,3},number_nodes::Int64,number_snapshots::Int64,α::Float64,R::Float64,ϵ::Float64)::Tuple{Array{Float64, 3},Array{Float64, 3}}

Update a set of nodes given a velocity eps.
"""
function updateNodes(coordinates::Array{Float64, 3},probabilities::Array{Float64,2},number_nodes::Int64,number_snapshots::Int64,α::Float64,R::Float64,ϵ::Float64)::Tuple{Array{Float64, 3},Array{Float64, 2}}
for t in 2:number_snapshots
    for i in 1:number_nodes
        p=probabilities[i,t-1]+rand(Uniform(-ϵ,ϵ))
        if p>1
            probabilities[i,t]=1-p%1
        elseif p<0
            probabilities[i,t]=abs(p)
        else
            probabilities[i,t]=p
        end
        coordinates[i,1,t]=(1/α)*acosh(1+(cosh(α*R)-1)*probabilities[i,t])
        coordinates[i,2,t]=coordinates[i,2,t-1]+rand(Uniform(-ϵ,ϵ))
    end
end
    coordinates, probabilities
end

function fillAdjacencyMatrix(adjacency::Array{Int64, 3},coordinates::Array{Float64, 3},number_nodes::Int64,number_snapshots::Int64,R::Float64,ζ::Float64)::Array{Int64, 3}
    for t in 1:number_snapshots
        for j in 1:number_nodes
            for i in 1:number_nodes
                if distanceHY(coordinates[i,:,t],coordinates[j,:,t],ζ)<=R
                    adjacency[i,j,t]=1
                end
            end
        end
    end
    return adjacency
end




#bad code/test 
function distanceHYTaylor(nodeA::Array{Float64, 1},nodeB::Array{Float64, 1}, ζ=1)::Float64
    hs(x)=x+(x^3)/factorial(3)+(x^5)/factorial(5)+(x^7)/factorial(7)+(x^9)/factorial(9)+(x^11)/factorial(11)
    hc(x)=x+(x^2)/factorial(2)+(x^4)/factorial(4)+(x^6)/factorial(6)+(x^8)/factorial(8)+(x^10)/factorial(10)
    println(cosh(ζ*nodeA[1]))
    println(hc(ζ*nodeA[1]))
    println(hs(ζ*nodeA[1]))
    if nodeA!=nodeB
        d=acosh(hc(ζ*nodeA[1])*hc(ζ*nodeB[1])-(hs(ζ*nodeA[1])*hs(ζ*nodeB[1])*cos(pi-abs(pi-abs(nodeA[2]-nodeB[2])))))/ζ
    else
        d=0.0
    end
    return d
end
function initializeNodesGrid(coordinates::Array{Float64, 3},probabilities::Array{Float64,2},number_nodes::Int64,α::Float64,R::Float64)::Tuple{Array{Float64, 3},Array{Float64, 2}}
    probabilities[:,1]=append!(collect(0.1*ones(number_nodes÷5)),collect(0.3*ones(number_nodes÷5)),collect(0.5*ones(number_nodes÷5)),collect(0.7*ones(number_nodes÷5)),collect(0.99*ones(number_nodes÷5)))
    coordinates[:,2,1]=append!(collect(LinRange(0,2*pi,number_nodes÷5)),collect(LinRange(0,2*pi,number_nodes÷5)),collect(LinRange(0,2*pi,number_nodes÷5)),collect(LinRange(0,2*pi,number_nodes÷5)),collect(LinRange(0,2*pi,number_nodes÷5)))
    for i in 1:number_nodes
        coordinates[i,1,1]=(1/α)*acosh(1+(cosh(α*R)-1)*probabilities[i,1])
        
    end
    return coordinates, probabilities
end

function initializeNodesGrid2(coordinates::Array{Float64, 3},probabilities::Array{Float64,2},number_nodes::Int64,α::Float64,R::Float64)::Tuple{Array{Float64, 3},Array{Float64, 2}}
    grid=collect(Iterators.product(1:sqrt(number_nodes), 1:sqrt(number_nodes)))
    k=1
    for j in 1:trunc(Int,sqrt(number_nodes))
        for i in 1:trunc(Int,sqrt(number_nodes))
            point=grid[i,j]
            probabilities[k,1]=rand(Uniform(0,1))
            coordinates[k,1,1]=sqrt(point[1]^2+point[2]^2)
            coordinates[k,2,1]=atan(point[2]/point[1])
            k+=1
            end
    end
    return coordinates, probabilities
end

struct HyperbolicTemporalNetworkTestGrid <: TemporalNetwork
    coordinates::Array{Float64, 3}
    probabilities::Array{Float64, 2}
    adjacency::Array{Int64, 3}
    number_nodes::Int64
    number_snapshots::Int64
    α::Float64
    R::Float64
    ϵ::Float64
    ζ::Float64

    function HyperbolicTemporalNetworkTestGrid(number_nodes::Int64,number_snapshots::Int64,α::Float64,R::Float64,ϵ::Float64,ζ::Float64)
        coordinates=Array{Float64}(undef, number_nodes,2,number_snapshots)
        probabilities=Array{Float64}(undef, number_nodes,number_snapshots)
        adjacency=fill(0,(number_nodes,number_nodes,number_snapshots))

        coordinates,probabilities=initializeNodesGrid(coordinates,probabilities,number_nodes,α,R)
        coordinates,probabilities=updateNodes(coordinates,probabilities,number_nodes,number_snapshots,α,R,ϵ)
        adjacency=fillAdjacencyMatrix(adjacency,coordinates,number_nodes,number_snapshots,R,ζ)
        new(coordinates,probabilities,adjacency,number_nodes,number_snapshots,α,R,ϵ,ζ)
    end
end