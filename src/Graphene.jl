module Graphene

using CSV: read
using Images: load, Gray, erode, label_components, component_centroids
using HDF5: h5read
using DataFrames
using NearestNeighbors: KDTree, knn, inrange
using GeometricalPredicates: Line2D

export CAtom, Atom2D, getx, gety, getxy, geti, gete

abstract type AbstractPoint2D end
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
getxy(g::Atom2D) = (getx(g), gety(g))
geti(p::Atom2D) = p._index
gete(p::Atom2D) = p._element

Bond(a::T, b::T) where {T<:Atom2D} = Line2D(a, b)
Bond(a::T, b::T) where {T<:AbstractPoint2D} = Line2D(a, b)
getx(l::Line2D) = (getx(l._a)+getx(l._b))/2
gety(l::Line2D) = (gety(l._a)+gety(l._b))/2
getxy(g::Line2D) = (getx(g), gety(g))
getl(l::Line2D) = (l._bx^2 + l._by^2)^0.5

include("fileio.jl")
include("atoms.jl")
include("bonds.jl")
include("defects.jl")


end # module
