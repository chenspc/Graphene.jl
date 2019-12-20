using Graphene
using Test

@test isequal([Tuple(test_xy[:,1])], index2xy(1, test_indexed_atoms))
@test isequal([Tuple(test_xy[:,12])], index2xy(12, test_indexed_atoms))
@test isequal(map(x -> Tuple(x), eachcol(test_xy[:,1:2:9])), index2xy(1:2:9, test_indexed_atoms))
