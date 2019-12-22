module Graphene

using CSV: read
using Images: load, Gray, erode, label_components, component_centroids
using HDF5: h5read
using DataFrames
using NearestNeighbors: KDTree, knn, inrange
using GeometricalPredicates: AbstractPoint2D, Point2D,Line, orientation
using JuliaDB: IndexedTable
import GeometricalPredicates: getx, gety, geta, getb

export CAtom, Atom2D, getx, gety, getxy, geti, gete, getl
export GAtom, GBond, GPolygon, GDefect

# struct Atom2D <: AbstractPoint2D
#     _x      ::Float64
#     _y      ::Float64
#     # _index  ::Int64
#     # _element::String
# end

# CAtom(x::Real, y::Real, index::Integer) = Atom2D(x, y, index, "C")
# CAtom() = CAtom(0., 0., 0)
# getx(p ::Atom2D) = p._x
# gety(p ::Atom2D) = p._y
# getxy(p::Atom2D) = (getx(p), gety(p))
# geti(p ::Atom2D) = p._index
# gete(p ::Atom2D) = p._element

# Bond(a::T, b::T) where {T<:Atom2D} = Line2D(a, b)
# getx(l ::Line2D{Atom2D}) = (getx(l._a)+getx(l._b))/2
# gety(l ::Line2D{Atom2D}) = (gety(l._a)+gety(l._b))/2
# getxy(l::Line2D{Atom2D}) = (getx(l), gety(l))
# geti(l ::Line2D{Atom2D}) = (geti(l._a), geti(l._b))
# gete(l ::Line2D{Atom2D}) = (gete(l._a), gete(l._b))
# getl(l ::Line2D{Atom2D}) = (l._bx^2 + l._by^2)^0.5

# getx(l) = (getx(l._a)+getx(l._b))/2
# gety(l) = (gety(l._a)+gety(l._b))/2
# getxy(l) = (getx(l), gety(l))
# geti(l) = (geti(l._a), geti(l._b))
# gete(l) = (gete(l._a), gete(l._b))
# getl(l) = (l._bx^2 + l._by^2)^0.5

abstract type AbstractGEntry end
mutable struct GAtom <: AbstractGEntry
    _id       ::UInt32
    _x        ::Float64
    _y        ::Float64
    _relatives::Vector{UInt32}
    _signature::String
    _frame    ::UInt64
    _dataset  ::String
    #
    # GAtom(id::UInt32, x::Float64, y::Float64, relatives::String, signature::String, frame::UInt64, dataset::String) = new(id, x, y, relatives, signature, frame, dataset)
end

mutable struct GBond <: AbstractGEntry
    _id       ::UInt32
    _x        ::Float64
    _y        ::Float64
    _relatives::Vector{UInt32}
    _signature::String
    _frame    ::UInt64
    _dataset  ::String
    #
    # GBond(id::UInt32, x::Float64, y::Float64, relatives::String, signature::String, frame::UInt64, dataset::String) = new(id, x, y, relatives, signature, frame, dataset)
end

mutable struct GPolygon <: AbstractGEntry
    _id       ::UInt32
    _x        ::Float64
    _y        ::Float64
    _relatives::Vector{UInt32}
    _signature::String
    _frame    ::UInt64
    _dataset  ::String

    _noa      ::UInt32
    #
    # GPolygon(id::UInt32, x::Float64, y::Float64, relatives::String, signature::String, frame::UInt64, dataset::String, noa::UInt32) = new(id, x, y, relatives, signature, frame, dataset, noa)
end

mutable struct GDefect<: AbstractGEntry
    _id       ::UInt32
    _x        ::Float64
    _y        ::Float64
    _relatives::Vector{UInt32}
    _signature::String
    _frame    ::UInt64
    _dataset  ::String

    _type     ::String
    #
    # GDefect(id::UInt32, x::Float64, y::Float64, relatives::String, signature::String, frame::UInt64, dataset::String, type::String) = new(id, x, y, relatives, signature, frame, dataset, type)
end

include("fileio.jl")
include("atoms.jl")
include("bonds.jl")
include("polygons.jl")
include("defects.jl")

end # module
