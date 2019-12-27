using Graphene
using Test

export test_bonds

@test isequal(path2turn((10, 11, 12)), (10, 11) => (11, 12))

test_bonds, test_paths = collect_bonds_paths(test_atom_group_collection, test_indexed_atoms)
test_turn_dict = Dict(map(path2turn, collect(test_paths)))
initial_length = length(test_turn_dict)
test_patoms, test_pbonds = link_turns!(test_turn_dict)
@test length(test_patoms) == length(test_pbonds) == initial_length - length(test_turn_dict)


test_polygons = make_polygons(test_paths)
@test length(test_polygons) == 3
@test Set(test_polygons[1][1]) == Set([5, 7, 9, 10, 11, 12, 13, 12, 14, 12, 11, 10, 8, 6, 4, 2, 1, 3])
@test Set(test_polygons[2][1]) == Set([7, 4, 6, 8, 10, 9])
@test Set(test_polygons[3][1]) == Set([1, 2, 4, 7, 5, 3])
