module Graphene

using StatsBase: countmap
using CSV: read
using Images: load, Gray, erode, label_components, component_centroids
using HDF5: h5read
using DataFrames
using NearestNeighbors: KDTree, knn, inrange
using GeometricalPredicates: Point2D, Line, orientation, getx, gety
using JuliaDB: IndexedTable
using AutoHashEquals
using Statistics: mean
using TupleTools
using SparseArrays: sparse

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

abstract type AbstractGPrimitive <: AbstractGEntry end
abstract type AbstractGDefect <: AbstractGEntry end

@auto_hash_equals mutable struct GAtom <: AbstractGPrimitive
    _id       ::Int
    _x        ::Float64
    _y        ::Float64
    _relatives::Set{Int}
    _signature::String
    _frame    ::Int
    _dataset  ::String
end
GAtom(id::Int, x::Float64, y::Float64, frame::Int, dataset::String) = GAtom(id, x, y, Set([]), "", frame, dataset)
GAtom(id::Int, x::Float64, y::Float64, frame::Int) = GAtom(id, x, y, Set([]), "", frame, "dataset")
GAtom(id::Int, x::Float64, y::Float64, relatives::Set) = GAtom(id, x, y, relatives, "", 0, "dataset")
GAtom(id::Int, x::Float64, y::Float64) = GAtom(id, x, y, Set([]), "", 0, "dataset")
GAtom(id::Int, p::Point2D) = GAtom(id, getx(p), gety(p), Set([]), "", 0, "dataset")
get_noa(g::GAtom) = 1
get_type(g::GAtom) = "Atom"

@auto_hash_equals mutable struct GBond <: AbstractGPrimitive
    _id       ::Int
    _x        ::Float64
    _y        ::Float64
    _relatives::Set{Int}
    _signature::String
    _frame    ::Int
    _dataset  ::String
end
GBond(id::Int, x::Float64, y::Float64, frame::Int, dataset::String) = GBond(id, x, y, Set([]), "", frame, dataset)
GBond(id::Int, x::Float64, y::Float64, frame::Int) = GBond(id, x, y, Set([]), "", frame, "dataset")
GBond(id::Int, x::Float64, y::Float64, relatives::Set) = GBond(id, x, y, relatives, "", 0, "dataset")
GBond(id::Int, x::Float64, y::Float64) = GBond(id, x, y, Set([]), "", 0, "dataset")
get_noa(g::GBond) = 2
get_type(g::GBond) = "Bond"

@auto_hash_equals mutable struct GPolygon <: AbstractGPrimitive
    _id       ::Int
    _x        ::Float64
    _y        ::Float64
    _relatives::Set{Int}
    _signature::String
    _frame    ::Int
    _dataset  ::String

    _noa      ::Int
end
GPolygon(id::Int, x::Float64, y::Float64, frame::Int, dataset::String, noa::Int) = GPolygon(id, x, y, Set([]), "", frame, dataset, noa)
GPolygon(id::Int, x::Float64, y::Float64, frame::Int, noa::Int) = GPolygon(id, x, y, Set([]), "", frame, "dataset", noa)
GPolygon(id::Int, x::Float64, y::Float64, noa::Int) = GPolygon(id, x, y, Set([]), "", 0, "dataset", noa)
get_noa(g::GPolygon) = g._noa
get_type(g::GPolygon) = "Polygon"

@auto_hash_equals mutable struct GDefect<: AbstractGDefect
    _id       ::Int
    _x        ::Float64
    _y        ::Float64
    _relatives::Set{Int}
    _signature::String
    _frame    ::Int
    _dataset  ::String

    _noa      ::Int
    _type     ::String
end
get_noa(g::GDefect) = g._noa
get_type(g::GDefect) = g._type

include("fileio.jl")
include("atoms.jl")
include("bonds.jl")
include("polygons.jl")
include("defects.jl")

end # module
