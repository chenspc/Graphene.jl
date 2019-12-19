# export Bond
# export Path
# export Polygon
# export my_isnothing
# export delete_nothing
# export atom_inrange
export make_atoms
export collect_atom_groups

function make_atoms(atom_xy)
    atoms_collection = map(CAtom, atom_xy[1,:], atom_xy[2,:], collect(1:size(atom_xy, 2)))
    indexed_atoms_collection = pairs(IndexLinear(), atoms_collection)
    return indexed_atoms_collection
end

function collect_atom_groups(atom_xy::Array{Float64,2}; max_bondlength=12)
    kdtree = KDTree(atom_xy, leafsize = 10)
    natoms = size(atom_xy, 2)
    indices, bondlengths = knn(kdtree, atom_xy, min(natoms, 4), true)
    indices_inrange = inrange(kdtree, atom_xy, max_bondlength, true)
    atoms_knn = zip(indices, indices_inrange)
    map(atom_inrange, atoms_knn)
end

function collect_atom_groups(atom_xy::Array{Float64,1}; max_bondlength=12)
    return [[1]]
end
#
function atom_inrange(atom_knn)
    index, index_inrange = atom_knn
    atom_selected = sort!(delete_nothing(indexin(index_inrange, index)))
    index_new = index[atom_selected]
    return index_new
end

function delete_nothing(x)
    if typeof(x) <: Array{Union{Nothing, T}} where T <: Number
        convert(Array{Int64}, filter!(my_isnothing, x))
    else
        x
    end
end

function my_isnothing(x)
    x != nothing
end
#
# # function graphene_model(args)
# #     body
# # end
