using Graphene
using GeometricalPredicates: Point2D, Line2D, Line
using Test

export test_xy
export test_indexed_atoms
export test_atom_group_collection

atom1 = Point2D(4., 6.)
atom2 = Point2D(1., 2.)
atom3 = Point2D(7., 2.)
atom4 = Point2D(8., 3.)
bond12 = Line(atom1, atom2)
@test isa(bond12, Line2D)

test_indexed_atoms = pairs(IndexLinear(), [atom1, atom2, atom3, atom4])
test_atom_group = [1, 2, 3, 4]
@test make_bonds(test_atom_group) == [(1,2), (1,3), (1,4)]
test_atom_group = [1, 2, 3]
@test make_bonds(test_atom_group) == [(1,2), (1,3)]
test_atom_group = [1, 2]
@test make_bonds(test_atom_group) == [(1,2)]

test_xy = [ 0.   5.    5.  15.   15.  20.  20.  30.  30.  35.  45.   55.   55.  64.   100.;
           17.4  8.7  26.1  8.7  26.1  0.  17.4  0.  17.4  8.7  8.7   8.7  18.7 10.4  100.]
test_indexed_atoms = make_atoms(test_xy)

test_atom_group = [4, 2, 6, 7]
@test make_paths(test_atom_group, test_indexed_atoms) == [(7, 4, 6), (6, 4, 2), (2, 4, 7)]
test_atom_group = [12, 11, 13, 14]
@test make_paths(test_atom_group, test_indexed_atoms) == [(11, 12, 13), (13, 12, 14), (14, 12, 11)]
test_atom_group = [11, 10, 12]
@test bond2path(test_atom_group, test_indexed_atoms) == ([(11, 10), (11, 12)], [(10, 11, 12), (12, 11, 10)])
test_atom_group = [15]
# @test bond2path(test_atom_group, test_indexed_atoms) == (nothing, nothing)
@test bond2path(test_atom_group, test_indexed_atoms) == ([], [])

test_atom_group_collection = [[1, 2, 3],
                              [2, 1, 4],
                              [3, 1, 5],
                              [4, 2, 6, 7],
                              [5, 3, 7],
                              [6, 4, 8],
                              [7, 4, 5, 9],
                              [8, 6, 10],
                              [9, 7, 10],
                              [10, 8, 9, 11],
                              [11, 10, 12],
                              [12, 11, 13, 14],
                              [13, 12],
                              [14, 12]]

@test collect_bonds_paths(test_atom_group_collection[9:14], test_indexed_atoms) ==
    ([(9, 7), (9, 10), (10, 8), (10, 9), (10, 11), (11, 10), (11, 12), (12, 11),
      (12, 13), (12, 14), (13, 12), (14, 12)],
     [(7, 9, 10), (10, 9, 7), (11, 10, 8), (8, 10, 9), (9, 10, 11),
     (10, 11, 12), (12, 11, 10), (11, 12, 13), (13, 12, 14), (14, 12, 11),
     (12, 13, 12), (12, 14, 12)])
