module Graphene

using CSV: read
using Images: load, Gray, erode, label_components, component_centroids
using HDF5: h5read
using DataFrames
using NearestNeighbors: KDTree, knn, inrange
using GeometricalPredicates: Point2D, Line, orientation, getx, gety
using JuliaDB: IndexedTable

# export getx, gety, getxy, geti, gete, getl
export GAtom, GBond, GPolygon, GDefect
export get_id, get_x, get_y, get_relatives, get_signature, get_frame, get_dataset, get_noa, get_type

abstract type AbstractGEntry end

get_id(g::AbstractGEntry) = g._id
get_x(g::AbstractGEntry) = g._x
get_y(g::AbstractGEntry) = g._y
get_relatives(g::AbstractGEntry) = g._relatives
get_signature(g::AbstractGEntry) = g._signature
get_frame(g::AbstractGEntry) = g._frame
get_dataset(g::AbstractGEntry) = g._dataset

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

get_noa(g::GAtom) = 1
get_type(g::GAtom) = "Atom"

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

get_noa(g::GBond) = 2
get_type(g::GBond) = "Bond"

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

get_noa(g::GPolygon) = g._noa
get_type(g::GPolygon) = "Polygon"

mutable struct GDefect<: AbstractGEntry
    _id       ::UInt32
    _x        ::Float64
    _y        ::Float64
    _relatives::Vector{UInt32}
    _signature::String
    _frame    ::UInt64
    _dataset  ::String

    _noa      ::UInt32
    _type     ::String
    #
    # GDefect(id::UInt32, x::Float64, y::Float64, relatives::String, signature::String, frame::UInt64, dataset::String, type::String) = new(id, x, y, relatives, signature, frame, dataset, type)
end

get_noa(g::GDefect) = g._noa
get_type(g::GDefect) = g._type

include("fileio.jl")
include("atoms.jl")
include("bonds.jl")
include("polygons.jl")
include("defects.jl")

end # module
