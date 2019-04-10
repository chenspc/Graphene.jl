using NearestNeighbors
using Test
using BenchmarkTools

@test my_isnothing(1) == true
@test my_isnothing(nothing) == false

a = [1, 3, 4, 5, 8, 11]
b = reverse!(collect(1:10))
ainb = indexin(a, b)

@test delete_nothing(a) == a
@test delete_nothing(ainb) == [10, 8, 7, 6, 3]

image_sampling = 1;
max_bondlength = 12.0;

kdtree = KDTree(atom_xy, leafsize = 10)
indices, bondlengths = knn(kdtree, atom_xy, 4, true)
indices_inrange = inrange(kdtree, atom_xy, max_bondlength, true)

f = zip(indices, bondlengths)
g = zip(indices, bondlengths, indices_inrange)
first(f)
first(g)

atom_inrange(first(g))
@time map(atom_inrange, g)
@benchmark map(atom_inrange, g)

@benchmark atom_groups(atom_xy)
atoms_info = atom_groups(atom_xy)
atom_info = atoms_info[4]
atoms, distances = atom_info
typeof(atoms_info)
eltype(atoms_info)
typeof(atom_info)
typeof(atoms)
typeof(distances)
length(atoms)

atoms0 = [1]
atoms1 = [1, 13]
atoms2 = [4, 9, 21]
atoms3 = [2, 8, 5, 22]

atoms3_xy = atom_xy[:,atoms3]
a0, a1, a2, a3 = map(Point, atoms3_xy[1,:], atoms3_xy[2,:])

@benchmark [(atoms3[1], atoms3[2]), (atoms3[1], atoms3[3]), (atoms3[1], atoms3[4])]
@benchmark deleteat!(collect(Iterators.product(first(atoms3), atoms3)), 1)

index, bondlength, index_inrange = first(g)

length(f)
typeof(f)
indices[1]
bondlengths[1]

indices_inrange[end]

c = indexin(indices_inrange[1], indices[1])
d = indexin(indices_inrange[end], indices[end])
e = indexin(indices_inrange[end-1], indices[end-1])

c_selected = sort!(delete_nothing(c))
d_selected = sort!(delete_nothing(d))
e_selected = sort!(delete_nothing(e))

indices[1][c_selected]
indices[end][d_selected]
indices[end-1][e_selected]

bondlengths[1][c_selected]
bondlengths[end][d_selected]
bondlengths[end-1][e_selected]


bond = Bond(CAtom(3.3,4.4,1), CAtom(5.5, 6.6,2))
orientation(bond, CAtom(1.2, 0.4, 3))
@benchmark length2(bond)

selected_atoms = atom_xy[:,atoms3]
a0, a1, a2, a3 = map(CAtom, selected_atoms[1,:], selected_atoms[2,:], atoms3)
b1, b2, b3 = broadcast(x -> Bond(a0, x), [a1, a2, a3])
map(orientation, [b1, b1, b2, b2, b3, b3], [a2, a3, a1, a3, a1, a2])
o12, o13, o21, o23, o31, o32 = map(orientation, [b1, b1, b2, b2, b3, b3], [a2, a3, a1, a3, a1, a2])


atoms3
atom_paths(atoms3, atom_xy)
@benchmark atom_paths(atoms3, atom_xy)

atoms3_a = [2, 5, 8, 22]
atoms3_b = [2, 22, 8, 5]
atoms3_c = [2, 22, 5, 8]
atom_paths(atoms3_a, atom_xy)
atom_paths(atoms3_b, atom_xy)
atom_paths(atoms3_c, atom_xy)
paths_a = Set(atom_paths(atoms3_a, atom_xy))
paths_b =Set(atom_paths(atoms3_b, atom_xy))
paths_c =Set(atom_paths(atoms3_c, atom_xy))
@test paths_a == paths_b && paths_a == paths_c && paths_b == paths_c
