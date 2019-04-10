# File test/runtests.jl
using Graphene
using Test

function tests()
    @testset "Subset of tests" begin
        include("test_id_generator.jl")
        include("test_data_inport.jl")
        # include("test_atoms.jl")
        # include("test_bonds.jl")
        # include("test_polygons.jl")
  end
end

@time tests()
