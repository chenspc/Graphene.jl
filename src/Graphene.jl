module Graphene

using StatsBase: counts
using CSV: read
using Images: load, Gray, erode, label_components, component_centroids
using HDF5: h5read
using DataFrames
using NearestNeighbors: KDTree, knn, inrange
using GeometricalPredicates: Point2D, Line, orientation, getx, gety, length2
using AutoHashEquals
using Statistics: mean
using SparseArrays: sparse, findnz
using Base.Iterators: Stateful, take, cycle
using Plots: plot, plot!, scatter, scatter!, @recipe, @colorant_str, @animate, gif
using Plots
using LazySets: convex_hull, VPolygon
using IterTools: subsets
using StructArrays

import Base.isless

export AbstractGEntry, AbstractGPrimitive, AbstractGDefect
export GAtom, GBond, GPolygon, GDefect
export get_id, get_cid, get_x, get_y, get_relatives, get_signature, get_frame, get_dataset, get_past, get_future
export get_noa, get_type, get_members
export isatom, isbond, ispolygon, istype, isfirst, islast, issingular

abstract type AbstractGEntry end
get_id(g::AbstractGEntry) = g._id
get_cid(g::AbstractGEntry) = g._id + g._frame * im
get_x(g::AbstractGEntry) = g._x
get_y(g::AbstractGEntry) = g._y
get_relatives(g::AbstractGEntry) = g._relatives
get_signature(g::AbstractGEntry) = g._signature
get_frame(g::AbstractGEntry) = g._frame
get_dataset(g::AbstractGEntry) = g._dataset
get_past(g::AbstractGEntry) = g._past
get_future(g::AbstractGEntry) = g._future

isatom(g::AbstractGEntry) = isequal(get_type(g), :atom)
isbond(g::AbstractGEntry) = isequal(get_type(g), :bond)
ispolygon(g::AbstractGEntry) = isequal(get_type(g), :polygon)
istype(g::AbstractGEntry, types...) = |([isequal(get_type(g), t) for t in types]...)
isless(a::AbstractGEntry, b::AbstractGEntry) = isless(get_noa(a), get_noa(b))
isfirst(a::AbstractGEntry) = isequal(get_past(a), 0im)
islast(a::AbstractGEntry) = isequal(get_future(a), 0im)
issingular(a::AbstractGEntry) = isfirst(a) && islast(a)
ispast(a::AbstractGEntry, b::AbstractGEntry) = isequal(get_cid(a), get_past(b))
isfuture(a::AbstractGEntry, b::AbstractGEntry) = isequal(get_cid(a), get_future(b))

abstract type AbstractGPrimitive <: AbstractGEntry end
get_members(g::AbstractGPrimitive) = Set(g._id)

abstract type AbstractGDefect <: AbstractGEntry end

@auto_hash_equals mutable struct GAtom <: AbstractGPrimitive
    _id       ::Int
    _x        ::Float64
    _y        ::Float64
    _relatives::Set{Int}
    _signature::String
    _frame    ::Int
    _dataset  ::String

    _past     ::Complex{Int}
    _future    ::Complex{Int}
end
GAtom(id::Int, x::Float64, y::Float64, relatives::Set, signature::String, frame::Int, dataset::String) = GAtom(id, x, y, relatives, signature, frame, dataset, 0im, 0im)
GAtom(id::Int, x::Float64, y::Float64, frame::Int, dataset::String) = GAtom(id, x, y, Set([]), "", frame, dataset)
GAtom(id::Int, x::Float64, y::Float64, frame::Int) = GAtom(id, x, y, frame, "dataset")
GAtom(id::Int, x::Float64, y::Float64, relatives::Set) = GAtom(id, x, y, relatives, "", 0, "dataset", 0im, 0im)
GAtom(id::Int, x::Float64, y::Float64) = GAtom(id, x, y, Set([]))
GAtom(id::Int, p::Point2D) = GAtom(id, getx(p), gety(p))
get_noa(g::GAtom) = 1
# get_type(g::GAtom) = :atom
get_type(g::GAtom) = :atom

@auto_hash_equals mutable struct GBond <: AbstractGPrimitive
    _id       ::Int
    _x        ::Float64
    _y        ::Float64
    _relatives::Set{Int}
    _signature::String
    _frame    ::Int
    _dataset  ::String

    _past     ::Complex{Int}
    _future    ::Complex{Int}
end
GBond(id::Int, x::Float64, y::Float64, relatives::Set, signature::String, frame::Int, dataset::String) = GBond(id, x, y, relatives, signature, frame, dataset, 0im, 0im)
GBond(id::Int, x::Float64, y::Float64, frame::Int, dataset::String) = GBond(id, x, y, Set([]), "", frame, dataset)
GBond(id::Int, x::Float64, y::Float64, frame::Int) = GBond(id, x, y, Set([]), "", frame, "dataset")
GBond(id::Int, x::Float64, y::Float64, relatives::Set) = GBond(id, x, y, relatives, "", 0, "dataset")
GBond(id::Int, x::Float64, y::Float64) = GBond(id, x, y, Set([]))
get_noa(g::GBond) = 2
# get_type(g::GBond) = :bond
get_type(g::GBond) = :bond

@auto_hash_equals mutable struct GPolygon <: AbstractGPrimitive
    _id       ::Int
    _x        ::Float64
    _y        ::Float64
    _relatives::Set{Int}
    _signature::String
    _frame    ::Int
    _dataset  ::String

    _noa      ::Int

    _past     ::Complex{Int}
    _future    ::Complex{Int}
end
GPolygon(id::Int, x::Float64, y::Float64, relatives::Set, signature::String, frame::Int, dataset::String, noa::Int) = GPolygon(id, x, y, relatives, signature, frame, dataset, noa, 0im, 0im)
GPolygon(id::Int, x::Float64, y::Float64, frame::Int, dataset::String, noa::Int) = GPolygon(id, x, y, Set([]), "", frame, dataset, noa)
GPolygon(id::Int, x::Float64, y::Float64, frame::Int, noa::Int) = GPolygon(id, x, y, frame, "dataset", noa)
GPolygon(id::Int, x::Float64, y::Float64, noa::Int) = GPolygon(id, x, y, 0, noa)
get_noa(g::GPolygon) = g._noa
# get_type(g::GPolygon) = :polygon
get_type(g::GPolygon) = :polygon

@auto_hash_equals mutable struct GEntry<: AbstractGEntry
    _id       ::Int
    _x        ::Float64
    _y        ::Float64
    _relatives::Set{Int}
    _signature::String
    _frame    ::Int
    _dataset  ::String

    _noa      ::Int
    _members  ::Set{Int}
    _type     ::Symbol

    _past     ::Complex{Int}
    _future    ::Complex{Int}
end
GEntry(id::Int64, x::Float64, y::Float64, relatives::Set, signature::String, frame::Int64, dataset::String, noa::Int64, members::Set) = GEntry(id, x, y, relatives, signature, frame, dataset, noa, members, :undetermined, 0im, 0im)
get_noa(g::GEntry) = g._noa
get_members(g::GEntry) = g._members
# get_type(g::GEntry) = "Undetermined"
# get_type(g::GEntry) = :undetermined

@auto_hash_equals mutable struct GDefect<: AbstractGDefect
    _id       ::Int
    _x        ::Float64
    _y        ::Float64
    _relatives::Set{Int}
    _signature::String
    _frame    ::Int
    _dataset  ::String

    _noa      ::Int
    _members  ::Set{Int}
    # _type     ::String
    _type     ::Symbol

    _past     ::Complex{Int}
    _future    ::Complex{Int}
end
GDefect(id::Int, x::Float64, y::Float64, relatives::Set, signature::String, frame::Int, dataset::String, noa::Int, members::Set{Int}, type::Symbol) = GDefect(id, x, y, relatives, signature, frame, dataset, noa, members, type, 0im, 0im)
get_noa(g::GDefect) = g._noa
get_members(g::GDefect) = g._members
get_type(g::GDefect) = g._type


# abstract type AbstractGraphene end
#
# @auto_hash_equals mutable struct Graphene <: AbstractGPrimitive
#     # _id       ::Int
#     # _x        ::Float64
#     # _y        ::Float64
#     # _relatives::Set{Int}
#     # _signature::String
#     _frame    ::Int
#     _dataset  ::String
#     _parts    ::AbstractArray{AbstractGPrimitive}
#     _defects  ::AbstractArray{AbstractGDefect}
#     #
#     # _past     ::Complex{Int}
#     # _future    ::Complex{Int}
# end
#
include("fileio.jl")
include("atoms.jl")
include("bonds.jl")
include("polygons.jl")
include("defects.jl")
include("plotting.jl")
include("series.jl")

end # module
