using Graphene
using Test
using BenchmarkTools

@test my_isnothing(1) == true
@test my_isnothing(nothing) == false

a = [1, 3, 4, 5, 8, 11]
b = reverse!(collect(1:10))
ainb = indexin(a, b)

@test delete_nothing(a) == a
@test delete_nothing(ainb) == [10, 8, 7, 6, 3]

atom_groups_collection = collect_atom_groups(atom_xy; max_bondlength=10.0)
# @benchmark atom_groups_collection = collect_atom_groups(atom_xy; max_bondlength=12.0)
atom_group = atom_groups_collection[4]
typeof(atom_groups_collection)
typeof(atom_group)
length(atom_group)

atoms0 = [1]
atoms1 = [1, 13]
atoms2 = [4, 9, 21]
atoms3 = [2, 8, 5, 22]

make_paths(atoms3, indexed_atoms_collection)
# @benchmark make_paths(atoms3, indexed_atoms_collection)
# @benchmark [(atoms3[1], atoms3[2]), (atoms3[1], atoms3[3]), (atoms3[1], atoms3[4])]
# @benchmark deleteat!(collect(Iterators.product(first(atoms3), atoms3)), 1)

bond = Bond(CAtom(3.3,4.4,1), CAtom(5.5, 6.6,2))

atoms3_a = [2, 5, 8, 22]
atoms3_b = [2, 5, 22, 8]
atoms3_c = [2, 8, 5, 22]
atoms3_d = [2, 8, 22, 5]
atoms3_e = [2, 22, 8, 5]
atoms3_f = [2, 22, 5, 8]
make_paths(atoms3_a, indexed_atoms_collection)
make_paths(atoms3_b, indexed_atoms_collection)
make_paths(atoms3_c, indexed_atoms_collection)
make_paths(atoms3_d, indexed_atoms_collection)
make_paths(atoms3_e, indexed_atoms_collection)
make_paths(atoms3_f, indexed_atoms_collection)
paths_a = Set(make_paths(atoms3_a, indexed_atoms_collection))
paths_b = Set(make_paths(atoms3_b, indexed_atoms_collection))
paths_c = Set(make_paths(atoms3_c, indexed_atoms_collection))
paths_d = Set(make_paths(atoms3_d, indexed_atoms_collection))
paths_e = Set(make_paths(atoms3_e, indexed_atoms_collection))
paths_f = Set(make_paths(atoms3_f, indexed_atoms_collection))
@test paths_a == paths_b == paths_c == paths_d == paths_e == paths_f

bonds_paths(atoms3_a, indexed_atoms_collection)

bondsNpaths = map(x -> bonds_paths(x, indexed_atoms_collection), atom_groups_collection)
# @benchmark bondsNpaths = map(x -> bonds_paths(x, indexed_atoms_collection), atom_groups_collection)
bonds = first.(bondsNpaths)
paths = last.(bondsNpaths)

bonds_collection, paths_collection = collect_bonds_paths(atom_groups_collection, indexed_atoms_collection)
# @benchmark bonds_collection, paths_collection = collect_bonds_paths(atom_groups_collection, indexed_atoms_collection)


function tally(zd, K)
    ret = zeros(Int64, K)
    for k in zd
        if k <= K
            ret[k] += 1
        end
    end
    return ret
end

atom_groups_collection = collect_atom_groups(atom_xy; max_bondlength=9.0)
bonds_collection, paths_collection = collect_bonds_paths(atom_groups_collection, indexed_atoms_collection)
polygon_collection = collect_polygons(paths_collection)
length.(first.(polygon_collection))
patoms_collection = first.(polygon_collection)
pbonds_collection = last.(polygon_collection)
# @benchmark atom_groups_collection = collect_atom_groups(atom_xy; max_bondlength=9.0)
# @benchmark onds_collection, paths_collection = collect_bonds_paths(atom_groups_collection, indexed_atoms_collection)
# @benchmark polygon_collection = collect_polygons(paths_collection)
# @benchmark tally(length.(first.(polygon_collection)), 20)
# @benchmark maximum(length.(first.(polygon_collection)))

tally(length.(first.(polygon_collection)), 20)
