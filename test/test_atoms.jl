using Graphene
using Test
using BenchmarkTools
using Statistics

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

atom_groups_collection = collect_atom_groups(atom_xy; max_bondlength=10.0)
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

# function markthebench(atom_xy)
#     atom_groups_collection = collect_atom_groups(atom_xy; max_bondlength=9.0)
#     bonds_collection, paths_collection = collect_bonds_paths(atom_groups_collection, indexed_atoms_collection)
#     polygon_collection = collect_polygons(paths_collection)
#     length.(first.(polygon_collection))
#     patoms_collection = first.(polygon_collection)
#     pbonds_collection = last.(polygon_collection)
# end
#
# markthebench(atom_xy)
# @benchmark markthebench(atom_xy)

patoms_collection
map(length, patoms_collection)
plocation_collection = map(length, patoms_collection)

plocation_collection = map(x->index2xy(x, indexed_atoms_collection), patoms_collection)

minimum(plocation_collection)

argmin(plocation_collection)

patoms_collection[48]

index2xy(patoms_collection[48],indexed_atoms_collection)
typeof(index2xy(patoms_collection[48],indexed_atoms_collection))
typeof(index2xy(patoms_collection[47],indexed_atoms_collection))
typeof(patoms_collection[48])
length(patoms_collection[48])

x = patoms_collection[48]

indexed_atoms_collection[patoms_collection[48]]

mean.([first.(plocation_collection[1]), last.(plocation_collection[1])])


plocation_collection[1]

patomsxy_collection = map(x->index2xy(x, indexed_atoms_collection), patoms_collection)
patomsxy_collection

px_collection = map(mean, map(x -> first.(x), patomsxy_collection))
py_collection = map(mean, map(x -> last.(x), patomsxy_collection))

zip(px_collection, py_collection)

px_collection[48]
py_collection[48]

patomsxy_collection[48]

pbonds_collection
pbonds_collection_sorted = pbonds_collection

pbonds1 = pbonds_collection[1]
pbonds2 = pbonds_collection[48]

collect(pbonds_collection)

collect(Iterators.flatten(pbonds_collection))

# a2b_dict, a2p_dict, b2a_dict, b2p_dict, a2bp_dict, b2ap_dict, p2ab_dict = make_relatives(adict, bdict, pdict, bond_dict)

patoms_collection
pnoa_collection = map(length, patoms_collection)
polygon_collection

Bond(indexed_atoms_collection[1], indexed_atoms_collection[2])

patoms_collection = first.(polygon_collection)
pbonds_collection = last.(polygon_collection)

# pnoa_collection = map(length, patoms_collection)
patomsxy_collection = map(x->index2xy(x, indexed_atoms_collection), patoms_collection)
px_collection = map(mean, map(x -> first.(x), patomsxy_collection))
py_collection = map(mean, map(x -> last.(x), patomsxy_collection))

all_bonds_directional = collect(Iterators.flatten(pbonds_collection))
f = x -> x=>Tuple(sort(collect(x)))
bond_dict = Dict(map(f, all_bonds_directional))
all_bonds = unique(collect(values(bond_dict)))
all_bonds1 = all_bonds[1]
testb = Bond(indexed_atoms_collection[first(all_bonds1)], indexed_atoms_collection[last(all_bonds1)])

getx(testb._a)
getx(testb)
gety(testb)
getxy(testb)

# pop!(unique!(collect(Iterators.flatten(filter!(x -> x != 0, map(x -> get(b2p_dict, x, 0), testv))))), testk)

g = make_graphene(data)

filter(x -> get(g._signature[x],7,0)==2 && get(g._signature[x],5,0)==2 && get(g._signature[x],6,0)==3, g._id)
filter(x -> g._noa[x]==12, g._id)
filter(x -> get(g._signature[x],7,0)==3, g._id)
filter(x -> get(g._signature[x],6,0)==6, g._id)
