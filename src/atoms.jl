# using DataFrames
# using DataFramesMeta
# using NearestNeighbors
# using GeometricalPredicates

# export CAtom, Atom2D, getx, gety, geti, gete
export Bond
export Path
export Polygon
export my_isnothing
export delete_nothing
export xy2atom
export index2xy
export atom_inrange
export collect_atom_groups
export make_bonds
export make_paths
export bonds_paths
export collect_bonds_paths
export make_path_dict
export make_polygon!
export collect_polygons


function my_isnothing(x)
    x != nothing
end

function delete_nothing(x)
    if typeof(x) <: Array{Union{Nothing, T}} where T <: Number
        convert(Array{Int64}, filter!(my_isnothing, x))
    else
        x
    end
end

function xy2atom(atom_xy)
     map(CAtom, atom_xy[1,:], atom_xy[2,:], collect(1:size(atom_xy, 2)))
end

function index2xy_old(x, indexed_atoms_collection)
    # if length(x) == 1
    #     atom = indexed_atoms_collection[x]
    #     output = getx(atom), gety(atom)
    # else
        atoms = indexed_atoms_collection[x]
        output = map(xx -> (getx(xx), gety(xx)), atoms)
    # end
    return output
end

function index2xy(x, indexed_atoms_collection)
    if length(x) == 1
        atom = indexed_atoms_collection[x[1]]
        output = [(getx(atom), gety(atom))]
    else
        atoms = indexed_atoms_collection[x]
        output = map(xx -> (getx(xx), gety(xx)), atoms)
    end
    return output
end

function atom_inrange(atom_knn)
    # index, bondlength, index_inrange = atom_knn
    index, index_inrange = atom_knn
    atom_selected = sort!(delete_nothing(indexin(index_inrange, index)))
    index_new = index[atom_selected]
    # bondlength_new = bondlength[atom_selected]
    # return index_new, bondlength_new
    return index_new
end

function collect_atom_groups(atom_xy; max_bondlength=12)
    kdtree = KDTree(atom_xy, leafsize = 10)
    indices, bondlengths = knn(kdtree, atom_xy, 4, true)
    indices_inrange = inrange(kdtree, atom_xy, max_bondlength, true)
    # atoms_knn = zip(indices, bondlengths, indices_inrange)
    atoms_knn = zip(indices, indices_inrange)
    map(atom_inrange, atoms_knn)
end

function make_bonds(atom_group)
    Set(deleteat!(collect(Iterators.product(first(atom_group), atom_group)), 1))
end

function make_paths(atoms3, indexed_atoms_collection)
    a0, a1, a2, a3 = indexed_atoms_collection[atoms3]
    b1, b2, b3 = broadcast(x -> Bond(a0, x), [a1, a2, a3])
    # o12, o13, o23, o21, o31, o32 = map(orientation, [b1, b1, b2, b2, b3, b3], [a2, a3, a3, a1, a1, a2])
    o12, o13, o23, o21, o31, o32 = map(GeometricalPredicates.orientation, [b1, b1, b2, b2, b3, b3], [a2, a3, a3, a1, a1, a2])
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
        nothing
    end

    # paths = map(x->Tuple(map(geti,x)), [[ps1 a0 pe1], [ps2 a0 pe2], [ps3 a0 pe3]])
    paths = Set(map(x->Tuple(map(geti,x)), [[ps1 a0 pe1], [ps2 a0 pe2], [ps3 a0 pe3]]))

end

function bonds_paths(atom_group::Array{Int64,1}, indexed_atoms_collection)
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

function collect_bonds_paths(atom_groups_collection, indexed_atoms_collection)
    bonds_paths_collection = map(x -> bonds_paths(x, indexed_atoms_collection), atom_groups_collection)
    bonds_collection = first.(bonds_paths_collection)
    paths_collection = last.(bonds_paths_collection)
    bonds_collection = union(filter!(x -> my_isnothing(x), bonds_collection)...)
    paths_collection = union(filter!(x -> my_isnothing(x), paths_collection)...)
    return bonds_collection, paths_collection
end

function make_path_dict(path)
    (path[1], path[2]) => (path[2], path[3])
end

function make_polygon!(path_dict)
    # carrier = first(first(path_dict))
    carrier = pop!(path_dict, first(first(path_dict)), "empty")
    patoms = [first(carrier)]
    pbonds = [carrier]
    while (carrier = pop!(path_dict, carrier, "empty")) != "empty"
        push!(patoms, first(carrier))
        push!(pbonds, carrier)
    end
    return patoms, pbonds
end

function collect_polygons(paths_collection)
    path_dict = Dict(map(make_path_dict, collect(paths_collection)))
    collector = []
    while !isempty(path_dict)
        push!(collector, make_polygon!(path_dict))
    end
    return collector
end
