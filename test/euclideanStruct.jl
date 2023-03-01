using Test, Statistics, Random
include("../src/euclideanModelStruct.jl")
include("../src/measuresStruct.jl") 

s=rand(1:1000)
Random.seed!(s)
println("seed ", s)

@testset verbose=true "Euclidean model" begin
    number_nodes=300
    number_snapshots=1
    d=140
    ϵ=0.06
    @testset verbose=true "Distance" begin
        @test distanceEU([1.0,1.0],[1.0,1.0])==0
        @test distanceEU([1.0,0.0],[0.0,1.0])==0
        @test distanceEU([0.0,0.0],[1.0,1.0])==0
        @test distanceEU([0.5,0.0],[0.5,1.0])==0
        @test distanceEU([0.23,0.11],[0.98,0.55])==distanceEU([0.98,0.55],[0.23,0.11])
    end
    @testset verbose=true "Degree" begin
 
        EUnetwork=EuclideanTemporalNetwork(number_nodes,number_snapshots,ϵ,d)
        @test averageDegree(EUnetwork)<((EUnetwork.r^2)*pi*(number_nodes-1)+2*sqrt((number_nodes-1)*pi*(EUnetwork.r^2)*(1-pi*(EUnetwork.r^2))))
        @test averageDegree(EUnetwork)>((EUnetwork.r^2)*pi*(number_nodes-1)-2*sqrt((number_nodes-1)*pi*(EUnetwork.r^2)*(1-pi*(EUnetwork.r^2))))
    end
    @testset verbose=true "UpdatingNodes" begin
        #test to check with Chebyshev inequality the radius 
        nodes=0.5*ones(100,2,2)
        nodes=updateNodes(nodes,100,2,ϵ)
        ρ=zeros(100)
        for i in 1:100
            ρ=sqrt(nodes[i,1,2]^2+nodes[i,2,2]^2)
            @test ρ<=ϵ+sqrt(0.5^2+0.5^2) #check that are in the circle
        end
        μ=mean(ρ)
        σ=var(ρ)
        λ=1.1
        outside=0
        inside=0
        for i in 1:100
            if ρ<μ-λ*σ || ρ>μ+λ*σ
                outside+=1
            else
                inside+=1
            end
        end
        @test outside/100<=(λ^2)
        @test inside/100>=1-(λ^2)

    end
end