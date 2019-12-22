using Graphene
using Test
using GeometricalPredicates:Point2D

# @test_broken CAtom() == Atom2D(0.0, 0.0, 0, "C")
# @test_broken test_atom == Atom2D(1.2, 1.5, 1, "C")

# @test_broken geti(test_atom) == 1
# @test_broken gete(test_atom) == "C"
# @test_broken getxy(test_atom) == (1.2, 1.5)

test_xy = [1. 2. 3.4;
           5  6  7.8]
test_atoms = make_atoms(test_xy)
# @test_broken test_atoms[1] == Atom2D(1.0, 5.0, 1, "C")
# @test_broken test_atoms[2] == Atom2D(2.0, 6.0, 2, "C")
# @test_broken test_atoms[3] == Atom2D(3.4, 7.8, 3, "C")

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

# delete_nothing()
