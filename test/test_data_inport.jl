using CSV
using DataFrames
using Test

export atom_xy
data = import_csv("/Users/chen/Dropbox/_julia/Graphene.jl/test/test_data.csv")
@test isa(data, DataFrames.DataFrame)
@test ncol(data) == 2

atom_xy = data_reshape(data)
@test size(atom_xy)[1] == 2
