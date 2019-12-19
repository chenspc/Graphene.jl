using Graphene
using Test

test_atom_groups = [[1]]
# bonds_paths(test_atom_groups)
bond2path(test_atom_groups)
test_atom_groups = [[1],
                    [2]]
test_atom_groups = [[1, 2],
                    [2, 1]]
test_atom_groups = [[1, 2, 3],
                    [2, 1, 4],
                    [3, 1, 5],
                    [4, 2, 6, 7],
                    [5, 3, 7],
                    [6, 4, 8],
                    [7, 4, 5, 9],
                    [8, 6, 10],
                    [9, 7, 10],
                    [10, 8, 9]]
