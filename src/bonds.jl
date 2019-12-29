export make_bonds
export make_paths
export bond2path
export collect_bonds_paths

function make_bonds(atom_group)
    Set(deleteat!(collect(Iterators.product(first(atom_group), atom_group)), 1))
end

function make_paths(atoms3, indexed_atoms)
    a0, a1, a2, a3 = indexed_atoms[atoms3]
    b1, b2, b3 = broadcast(x -> Line(a0, x), [a1, a2, a3])
    o12, o13, o23, o21, o31, o32 = map(orientation, [b1, b1, b2, b2, b3, b3], [a2, a3, a3, a1, a1, a2])
    p0, p1, p2, p3 = Tuple(atoms3)
    if     o12 * o13 == -1
        if o12 == -1
           ps1, pe1, ps2, pe2, ps3, pe3 = p3, p1, p1, p2, p2, p3
        else
           ps1, pe1, ps2, pe2, ps3, pe3 = p3, p2, p2, p1, p1, p3
        end
    elseif o23 * o21 == -1
        if o23 == -1
           ps1, pe1, ps2, pe2, ps3, pe3 = p1, p2, p2, p3, p3, p1
        else
           ps1, pe1, ps2, pe2, ps3, pe3 = p1, p3, p3, p2, p2, p1
        end
    elseif o31 * o32 == -1
        if o31 == -1
           ps1, pe1, ps2, pe2, ps3, pe3 = p2, p3, p3, p1, p1, p2
        else
           ps1, pe1, ps2, pe2, ps3, pe3 = p2, p1, p1, p3, p3, p2
        end
    else
        path_exist = 0
    end
    paths = Set(map(x->Tuple(x), [[ps1 p0 pe1], [ps2 p0 pe2], [ps3 p0 pe3]]))
end

function bond2path(atom_group::Array{Int64,1}, indexed_atoms)
    nn = length(atom_group) - 1
    if     nn == 3
        bonds = make_bonds(atom_group)
        paths = make_paths(atom_group, indexed_atoms)
    elseif nn == 2
        bonds = make_bonds(atom_group)
        path = (atom_group[2], atom_group[1], atom_group[3])
        paths = Set([path, reverse(path)])
    elseif nn == 1
        bonds = make_bonds(atom_group)
        paths = Set([(atom_group[2], atom_group[1], atom_group[2])])
    else   nn == 0
        bonds = Set([])
        paths = Set([])
    end
    return bonds, paths
end

function collect_bonds_paths(atoms_groups, indexed_atoms)
    bonds_paths_collection = map(x -> bond2path(x, indexed_atoms), atoms_groups)
    bonds = union(first.(bonds_paths_collection)...)
    paths = union(last.(bonds_paths_collection)...)
    return bonds, paths
end
