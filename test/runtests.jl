using Test

testsmodels=["euclideanStruct","hyperbolicStruct"]
tests=["measuresStruct"]
const testdir = dirname(@__FILE__)
@testset "All tests Struct" verbose=true begin
    @testset "Tools" verbose=true begin
        for t in tests
            tp = joinpath(testdir, "$(t).jl")
            include(tp)
        end
    end
    @testset "Graph models" verbose=true begin
        for t in testsmodels
            tp = joinpath(testdir, "$(t).jl")
            include(tp)
        end
    end
end