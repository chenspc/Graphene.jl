using Graphene
using Test
using GeometricalPredicates:Point2D

test_xy = [0;
           17.4]
test_group = Set.([[1]])
test_atom_groups = collect_atom_groups(test_xy)
@test test_group == Set.(test_atom_groups)

test_xy = [0     5;
           17.4  8.7]
test_group = Set.([[1, 2], [2, 1]])
test_atom_groups = collect_atom_groups(test_xy)
@test test_group == Set.(test_atom_groups)

test_xy = [0     35;
           17.4  8.7]
test_group = Set.([[1],[2]])
test_atom_groups = collect_atom_groups(test_xy)
@test test_group == Set.(test_atom_groups)

test_xy = [0     5     5   15    15   20   20   30   30   35;
           17.4  8.7  26.1  8.7  26.1  0   17.4  0   17.4  8.7]
test_group = Set.([[1, 2, 3],
                    [2, 1, 4],
                    [3, 1, 5],
                    [4, 2, 6, 7],
                    [5, 3, 7],
                    [6, 4, 8],
                    [7, 4, 5, 9],
                    [8, 6, 10],
                    [9, 7, 10],
                    [10, 8, 9]])
test_atom_groups = collect_atom_groups(test_xy)
@test test_group == Set.(test_atom_groups)
