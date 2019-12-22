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
