export make_bonds
export make_paths
export bond2path
export collect_bonds_paths

function make_bonds(atom_group)
    Set(deleteat!(collect(Iterators.product(first(atom_group), atom_group)), 1))
end

function make_paths(atoms3, indexed_atoms_collection)
    a0, a1, a2, a3 = indexed_atoms_collection[atoms3]
    b1, b2, b3 = broadcast(x -> Bond(a0, x), [a1, a2, a3])
    o12, o13, o23, o21, o31, o32 = map(orientation, [b1, b1, b2, b2, b3, b3], [a2, a3, a3, a1, a1, a2])
    # o12, o13, o23, o21, o31, o32 = map(GeometricalPredicates.orientation, [b1, b1, b2, b2, b3, b3], [a2, a3, a3, a1, a1, a2])
    path_exist = 1
    if     o12 * o13 == -1
        if o12 == -1
           ps1, pe1, ps2, pe2, ps3, pe3 = a3, a1, a1, a2, a2, a3
        else
           ps1, pe1, ps2, pe2, ps3, pe3 = a3, a2, a2, a1, a1, a3
        end
    elseif o23 * o21 == -1
        if o23 == -1
           ps1, pe1, ps2, pe2, ps3, pe3 = a1, a2, a2, a3, a3, a1
        else
           ps1, pe1, ps2, pe2, ps3, pe3 = a1, a3, a3, a2, a2, a1
        end
    elseif o31 * o32 == -1
        if o31 == -1
           ps1, pe1, ps2, pe2, ps3, pe3 = a2, a3, a3, a1, a1, a2
        else
           ps1, pe1, ps2, pe2, ps3, pe3 = a2, a1, a1, a3, a3, a2
        end
    else
        # nothing
        path_exist = 0
    end

    if path_exist == 1
        # paths = map(x->Tuple(map(geti,x)), [[ps1 a0 pe1], [ps2 a0 pe2], [ps3 a0 pe3]])
        paths = Set(map(x->Tuple(map(geti,x)), [[ps1 a0 pe1], [ps2 a0 pe2], [ps3 a0 pe3]]))
    else
        # paths = nothing
        paths = Set([])
    end

    return paths

end

function bond2path(atom_group::Array{Int64,1}, indexed_atoms_collection)
    # atoms, bond_lengths = atom_info
    nn = length(atom_group) - 1
    if     nn == 3
        bonds = make_bonds(atom_group)
        paths = make_paths(atom_group, indexed_atoms_collection)
    elseif nn == 2
        bonds = make_bonds(atom_group)
        path = (atom_group[2], atom_group[1], atom_group[3])
        paths = Set([path, reverse(path)])
    elseif nn == 1
        bonds = make_bonds(atom_group)
        paths = Set([(atom_group[2], atom_group[1], atom_group[2])])
    else   nn == 0
        bonds = nothing
        paths = nothing
    end
    return bonds, paths
end

function collect_bonds_paths(atom_group_collection, indexed_atoms_collection)
    bonds_paths_collection = map(x -> bond2path(x, indexed_atoms_collection), atom_group_collection)
    bonds_collection = first.(bonds_paths_collection)
    paths_collection = last.(bonds_paths_collection)
    bonds_collection = union(filter!(x -> my_isnothing(x), bonds_collection)...)
    paths_collection = union(filter!(x -> my_isnothing(x), paths_collection)...)
    return bonds_collection, paths_collection
end
