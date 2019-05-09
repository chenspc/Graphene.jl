using Graphene
using Test
using BenchmarkTools
using DataFrames

export atom_xy
export atoms_collection
export indexed_atoms_collection
# data = import_csv("/Users/chen/Dropbox/_julia/Graphene.jl/test/test_data.csv")
data = import_csv("/Users/chen/Dropbox/_julia/Graphene.jl/test/test_data_flower.csv")
@test isa(data, DataFrames.DataFrame)
@test ncol(data) == 2

atom_xy = data_reshape(data)
@test size(atom_xy)[1] == 2

# atoms_collection = xy2atom(atom_xy)
# indexed_atoms_collection = pairs(IndexLinear(), atoms_collection)
# @benchmark atoms_collection = xy2atom(atom_xy)
# @benchmark indexed_atoms_collection = pairs(IndexLinear(), atoms_collection)

indexed_atoms_collection = xy2atom(atom_xy)

# @test typeof(atoms_collection) == Array{Graphene.Atom2D,1}
# @test typeof(atoms_collection) == Array{Atom2D,1}
