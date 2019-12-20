module Graphene

using CSV: read
using Images: load, Gray, erode, label_components, component_centroids
using HDF5: h5read
using DataFrames
using NearestNeighbors: KDTree, knn, inrange
using GeometricalPredicates: AbstractPoint2D, Line2D, orientation
using JuliaDB: IndexedTable
import GeometricalPredicates: getx, gety, geta, getb

export CAtom, Atom2D, Bond, getx, gety, getxy, geti, gete, getl

struct Atom2D <: AbstractPoint2D
    _x::Float64
    _y::Float64
    _index::Int64
    _element::String
    Atom2D(x::Float64, y::Float64, index::Int64, element::String) = new(x, y, index, element)
end

CAtom(x::Real, y::Real, index::Integer) = Atom2D(x, y, index, "C")
CAtom() = CAtom(0., 0., 0)
getx(p::Atom2D) = p._x
gety(p::Atom2D) = p._y
getxy(p::Atom2D) = (getx(p), gety(p))
geti(p::Atom2D) = p._index
gete(p::Atom2D) = p._element

Bond(a::T, b::T) where {T<:Atom2D} = Line2D(a, b)
getx(l::Line2D{Atom2D}) = (getx(l._a)+getx(l._b))/2
gety(l::Line2D{Atom2D}) = (gety(l._a)+gety(l._b))/2
getxy(l::Line2D{Atom2D}) = (getx(l), gety(l))
geti(l::Line2D{Atom2D}) = (geti(l._a), geti(l._b))
gete(l::Line2D{Atom2D}) = (gete(l._a), gete(l._b))
getl(l::Line2D{Atom2D}) = (l._bx^2 + l._by^2)^0.5

include("fileio.jl")
include("atoms.jl")
include("bonds.jl")
include("polygons.jl")
include("defects.jl")


end # module
