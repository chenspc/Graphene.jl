module Graphene

using DataFrames
using DataFramesMeta
using JuliaDB
using OnlineStats
using StatsBase
using Query
using CSV
# using CategoricalArrays
using Missings
using IterTools
using GeometricalPredicates
using NearestNeighbors
using Images
# using ImageMorphology
# using ImageSegmentation
# using ImageView
using FileIO
using QuartzImageIO
using HDF5
# using Gadfly
using Plots
using Random
using Base64
# using Luxor
# using Colors
# using Images
using Statistics
using BenchmarkTools

export CAtom, Atom2D, getx, gety, geti, gete,
       Bond, getl,
       Path,
       # Polygon,
       AbstractGElement, GElement, getid, gettype, getxy, getrelatives, getsignature, getnoa, getframe, getdataset,
       AbstractGraphene, GFrame

import GeometricalPredicates: getx, gety, geta, getb

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
getxy(g::Atom2D) = (getx(g), gety(g))
geti(p::Atom2D) = p._index
gete(p::Atom2D) = p._element

CAtom(x::Real, y::Real, index::Integer) = Atom2D(x, y, index, "C")

Bond(a::T, b::T) where {T<:Atom2D} = Line2D(a, b)
Bond(a::T, b::T) where {T<:AbstractPoint2D} = Line2D(a, b)
# getx(l::Line2D) = l._bx/2
# gety(l::Line2D) = l._by/2
getx(l::Line2D) = (getx(l._a)+getx(l._b))/2
gety(l::Line2D) = (gety(l._a)+gety(l._b))/2
getxy(g::Line2D) = (getx(g), gety(g))
getl(l::Line2D) = (l._bx^2 + l._by^2)^0.5

abstract type AbstractGElement end

# struct GElement <: AbstractGElement
#     _id::UInt32
#     _noa::Int32
#     _x::Float64
#     _y::Float64
#     _relatives::Tuple{Vararg{UInt}}
#     # _signature::Tuple{Int, Tuple{Vararg{Int}}}
#     _signature::Dict{Int,Int}
#     _frame::Int64
#     _dataset::String
#     GElement(id::UInt32, noa::Int32, x::Float64, y::Float64, relatives::Tuple{Vararg{UInt}}, signature::Dict{Int,Int}, frame::Int64, dataset::String) = new(id, noa, x, y, relatives, signature, frame, dataset)
# end
# GElement() = GElement(hash(0), 0, 0., 0., (), Dict(), 0, "empty")
# getid(g::GElement) = g._id
# getnoa(g::GElement) = g._noa
# getx(g::GElement) = g._x
# gety(g::GElement) = g._y
# getxy(g::GElement) = (getx(g), gety(g))
# getrelatives(g::GElement) = g._relatives
# getsignature(g::GElement) = g._signature
# getframe(g::GElement) = g._frame
# getdataset(g::GElement) = g._dataset
#
abstract type AbstractGraphene end

struct GFrame <: AbstractGraphene
    _id::Array{UInt32,1}
    _type::Dict{UInt32,String}
    _x::Dict{UInt32,Float64}
    _y::Dict{UInt32,Float64}
    # _relatives::Dict{UInt32,Array{UInt32,1}}
    _relatives::Dict{UInt32,String}
    # _signature::Dict{UInt32, Dict{Int,Int}}
    _signature::Dict{UInt32, String}
    _noa::Dict{UInt32, Int64}
    _frame::Int64
    _dataset::String
    GFrame(id::Array{UInt32,1}, type::Dict{UInt32,String}, x::Dict{UInt32,Float64}, y::Dict{UInt32,Float64}, relatives::Dict{UInt32,String}, signature::Dict{UInt32, String}, noa::Dict{UInt32,Int64}, frame::Int64, dataset::String) = new(id, type, x, y, relatives, signature, noa, frame, dataset)
end
getid(g::GFrame) = g._id
gettype(g::GFrame) = g._type
getx(g::GFrame) = g._x
gety(g::GFrame) = g._y
getxy(g::GFrame) = (getx(g), gety(g))
getrelatives(g::GFrame) = g._relatives
getsignature(g::GFrame) = g._signature
getnoa(g::GFrame) = g._noa
getframe(g::GFrame) = g._frame
getdataset(g::GFrame) = g._dataset

# struct Bond2D{T<:AbstractPoint2D} <: Line2D
# struct Bond2D{T<:Atom2D} <: Line2D
#     _a::T
#     _b::T
#     _bx::Float64
#     _by::Float64
#     _x::Float64
#     _y::Float64
#     _length::Float64
#     function Bond2D{T}(a::T, b::T) where T
#         bx = getx(b) - getx(a)
#         by = gety(b) - gety(a)
#         x = bx/2
#         y = by/2
#         length = (bx^2 + by^2)^0.5
#         new(a, b, bx, by, x, y, length)
#     end
# end

# Bond2D(a::T, b::T) where {T<:Atom2D} = Bond2D{T}(a, b)



# Change this to enable debugging
const DEBUG = false

include("data_import.jl")
include("id_generator.jl")
include("points.jl")
include("atoms.jl")
include("plotting.jl")
include("defect.jl")

end # module
