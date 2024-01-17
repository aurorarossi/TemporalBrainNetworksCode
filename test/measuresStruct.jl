using Test, Random
include("../src/measuresStruct.jl")

s = rand(1:1000)
Random.seed!(s)
println("seed ", s)

abstract type TemporalNetwork end

struct TestTemporalNetwork <: TemporalNetwork
    adjacency::Array{Int64, 3}
    number_nodes::Int64
    number_snapshots::Int64
end

@testset verbose = true "Measures" begin 
    matrix = zeros(5, 5, 4)
    matrix[ 1, 2,1] = 1
    matrix[ 2, 1,1] = 1
    matrix[ 3, 4,1] = 1
    matrix[ 4, 3,1] = 1
    matrix[ 2, 3,2] = 1
    matrix[ 3, 2,2] = 1
    matrix[ 1, 5,3] = 1
    matrix[ 5, 1,3] = 1
    matrix[ 2, 5,4] = 1
    matrix[ 5, 2,4] = 1
    matrix[ 3, 5,4] = 1
    matrix[ 5, 3,4] = 1
    matrix[ 3, 4,4] = 1
    matrix[ 4, 3,4] = 1
    testnetwork=TestTemporalNetwork(matrix,5,4)
    BFSres = [[0.0, 1.0, 2.0, 4.0, 3.0], [1.0, 0.0, 2.0, 4.0, 3.0], [4.0, 2.0, 0.0, 1.0, 4.0], [4.0, 2.0, 1.0, 0.0, 4.0], [3.0, 4.0, 4.0, 4.0, 0.0]]
    matrixrand = deepcopy(matrix)
    for _ in 1:40
        i = rand(1:5)
        j = rand(1:5)
        t = rand(1:4)
        matrixrand[ i, j,t] = 1
        matrixrand[ j, i,t] = 1
    end
    
    testnetworkrand=TestTemporalNetwork(matrixrand,5,4)
    @testset verbose = true "BFS" begin
        for i in 1:5
            @test BFS(i, testnetwork) == BFSres[i]
        end
    end
    
    @testset verbose = true "temporalPathLength" begin
        @test temporalPathLength(testnetwork) == 57 / 20
        @test temporalPathLength(testnetwork) > temporalPathLength(testnetworkrand)
    end

    @testset verbose = true "temporalCorrelationCoefficient" begin
        matrix[ 3, 2,1] = 1
        matrix[ 2, 3,1] = 1
        testnetwork3=TestTemporalNetwork(matrix,5,4)
        @test temporalCorrelationCoefficient(testnetwork) == 0
        @test round(temporalCorrelationCoefficient(testnetwork3),digits=6) == round(2/(15*sqrt(2)),digits=6)
    end

    @testset verbose = true "clusteringCoeff" begin
        matrix2 = zeros(6, 6, 1)
        matrix2[ 1, 2,1] = 1
        matrix2[ 2, 1,1] = 1
        matrix2[ 1, 5,1] = 1
        matrix2[ 5, 1,1] = 1
        matrix2[ 2, 5,1] = 1
        matrix2[ 5, 2,1] = 1
        matrix2[ 2, 3,1] = 1
        matrix2[ 3, 2,1] = 1
        matrix2[ 4, 5,1] = 1
        matrix2[ 5, 4,1] = 1
        matrix2[ 4, 3,1] = 1
        matrix2[ 3, 4,1] = 1
        matrix2[ 4, 6,1] = 1
        matrix2[ 6, 4,1] = 1
        testnetwork2=TestTemporalNetwork(matrix2,6,1)
        @test clustringCoeffOneSnap(1, testnetwork) == 0
        @test clustringCoeffOneSnap(1, testnetwork2) == 3 / 11
        @test clustringCoeffOneSnap(1, testnetwork) <= clustringCoeffOneSnap(1, testnetworkrand)
        @test clustringCoeffOneSnap(1, testnetworkrand) >= 0 && clustringCoeffOneSnap(1, testnetworkrand) <= 1
    end

    @testset verbose = true "temporalClosenessCentrality" begin
        cCres=[(1/4)*(1+1/2+1/4+1/3),(1/4)*(1+1/2+1/4+1/3),(1/4)*(1/2+1+1/4),(1/4)*(1/2+1+1/4),(1/4)*(1/3+1/2)]
        for i in 1:4
            @test temporalClosenessCentrality(testnetwork,i,1)==cCres[i]
        end
    end
end


