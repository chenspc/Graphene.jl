using Graphene
using Test

function tests()
    @testset "Subset of tests" begin
        include("test_fileio.jl")
        include("test_atoms.jl")
  end
end

tests()
@time tests()
