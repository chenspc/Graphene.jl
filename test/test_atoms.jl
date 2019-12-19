using Graphene
using Test

@test CAtom() == Atom2D(0.0, 0.0, 0, "C")
test_atom = CAtom(1.2, 1.5, 1)
@test test_atom == Atom2D(1.2, 1.5, 1, "C")

@test getx(test_atom) == 1.2
@test gety(test_atom) == 1.5
@test geti(test_atom) == 1
@test gete(test_atom) == "C"
@test getxy(test_atom) == (1.2, 1.5)

test_xy = [1. 2. 3.4;
           5  6  7.8]
test_atoms = make_atoms(test_xy)
@test test_atoms[1] == Atom2D(1.0, 5.0, 1, "C")
@test test_atoms[2] == Atom2D(2.0, 6.0, 2, "C")
@test test_atoms[3] == Atom2D(3.4, 7.8, 3, "C")

test_xy = [0;
           17.3]
test_group = Set.([[1]])
test_atom_groups = collect_atom_groups(test_xy)
@test test_group == Set.(test_atom_groups)

test_xy = [0     5;
           17.3  8.7]
test_group = Set.([[1, 2], [2, 1]])
test_atom_groups = collect_atom_groups(test_xy)
@test test_group == Set.(test_atom_groups)

test_xy = [0     35;
           17.3  8.7]
test_group = Set.([[1],[2]])
test_atom_groups = collect_atom_groups(test_xy)
@test test_group == Set.(test_atom_groups)

test_xy = [0     5    5   15   15  20  20    30  30    35;
           17.3  8.7  26  8.7  26  0   17.3  0   17.3  8.7]
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
