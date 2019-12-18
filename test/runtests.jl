using Graphene
using Test

function tests()
    @testset "Subset of tests" begin
        include("test_fileio.jl")
  end
end

tests()
@time tests()
