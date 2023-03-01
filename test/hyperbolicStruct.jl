using Test, Random, Distributions
include("../src/hyperbolicModelStruct.jl")
include("../src/measuresStruct.jl")
s=rand(1:1000)
Random.seed!(s)
println("seed ", s)
@testset verbose = true "Hyperbolic model" begin
    number_nodes=50
    number_snapshots=23
    α = 2.0
    R = 1.0
    ζ = 1.0
    ϵ = 0.9

    @testset verbose=true "Distance" begin
        @test distanceHY([1.0,1.0],[1.0,1.0])==0
        @test distanceHY([0.23,0.11],[0.98,0.55])==distanceHY([0.98,0.55],[0.23,0.11])
    end



    @testset verbose=true "initializeNodes and updateNodes" begin
        testnetwork=HyperbolicTemporalNetwork(number_nodes,number_snapshots,α,R,ϵ,ζ)
        for i in 1:number_nodes
            @test testnetwork.probabilities[i, 1] >= 0 && testnetwork.probabilities[i, 1] <= 1
            @test testnetwork.probabilities[i,2,1] >= 0 && testnetwork.probabilities[i,2,1] <= 2 * pi
            @test testnetwork.probabilities[i,1,1] <= R

        end
    end

    @testset verbose=true "check updateNodes wrt speed" begin
        testnetwork1=HyperbolicTemporalNetwork(number_nodes,number_snapshots,α,R,0.01,ζ)
        testnetwork2=HyperbolicTemporalNetwork(number_nodes,number_snapshots,α,R,0.3,ζ)
        testnetwork3=HyperbolicTemporalNetwork(number_nodes,number_snapshots,α,R,0.8,ζ)
        @test temporalCorrelationCoefficient(testnetwork1)>=temporalCorrelationCoefficient(testnetwork2)
        @test temporalCorrelationCoefficient(testnetwork2)>=temporalCorrelationCoefficient(testnetwork3)


    end
end