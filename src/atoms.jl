using DataFrames
using DataFramesMeta
using NearestNeighbors
using GeometricalPredicates

export Atom, Atom2D, getx, gety, geti, gete
export Bond
export my_isnothing
export delete_nothing
export atom_inrange
export atom_groups

import GeometricalPredicates.getx
import GeometricalPredicates.gety

struct Atom2D <: AbstractPoint2D
    _x::Float64
    _y::Float64
    _index::Int64
    _element::String
    Atom2D(x::Float64, y::Float64, index::Int64, element::String) = new(x, y, index, element)
end
Atom2D() = Atom2D(0., 0., 0, "C")
getx(p::Atom2D) = p._x
gety(p::Atom2D) = p._y
geti(p::Atom2D) = p._index
gete(p::Atom2D) = p._element

CAtom(x::Real, y::Real, index::Integer) = Atom2D(x, y, index, "C")
Bond(a::T, b::T) where {T<:Atom2D} = Line2D(a, b)

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

function atom_inrange(atom_knn)
    index, bondlength, index_inrange = atom_knn
    atom_selected = sort!(delete_nothing(indexin(index_inrange, index)))
    index_new = index[atom_selected]
    bondlength_new = bondlength[atom_selected]
    return index_new, bondlength_new
end

function atom_groups(atom_xy; max_bondlength=12)
    kdtree = KDTree(atom_xy, leafsize = 10)
    indices, bondlengths = knn(kdtree, atom_xy, 4, true)
    indices_inrange = inrange(kdtree, atom_xy, max_bondlength, true)
    atoms_knn = zip(indices, bondlengths, indices_inrange)
    map(atom_inrange, atoms_knn)
end

function atom_paris(atoms)
    deleteat!(collect(Iterators.product(first(atoms), atoms)), 1)
end

function atom_paths(atoms, atom_xy)
    selected_atoms = atom_xy[:,atoms]
    a0, a1, a2, a3 = map(CAtom, selected_atoms[1,:], selected_atoms[2,:], atoms3)
    b1, b2, b3 = broadcast(x -> Bond(a0, x), [a1, a2, a3])
    o12, o13, o23, o21, o31, o32 = map(orientation, [b1, b1, b2, b2, b3, b3], [a2, a3, a3, a1, a1, a2])
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

    paths = map(x->Tuple(map(geti,x)), [[ps1 a0 pe1], [ps2 a0 pe2], [ps3 a0 pe3]])

end

function atoms2bonds(atom_info::Tuple{Array{Int64,1},Array{Float64,1}})
    atoms, bond_lengths = atom_info
    nn = length(atoms) - 1
    if     nn == 3
        bonds = Set(atom_paris(atoms))
        paths = atom_paths(atoms, atom_xy)
    elseif nn == 2
        bonds = Set(atom_paris(atoms))
        path = (atoms[2], atoms[1], atoms[3])
        paths = [path, reverse(path)]
    elseif nn == 1
        bonds = Set(atom_paris(atoms))
        paths = (atoms[2], atoms[1], atoms[2])
    else   nn == 0
        bonds = nothing
        paths = nothing
    end
end
