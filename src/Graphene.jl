module Graphene

using DataFrames
using DataFramesMeta
using Query
using CSV
using CategoricalArrays
using Missings
using IterTools
using GeometricalPredicates
using NearestNeighbors
# using Gadfly
using Plots
# using Luxor
# using Colors
# using Images
using BenchmarkTools

export CAtom, Atom2D, getx, gety, geti, gete,
       Bond, getl,
       Path,
       Polygon

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
geti(p::Atom2D) = p._index
gete(p::Atom2D) = p._element

CAtom(x::Real, y::Real, index::Integer) = Atom2D(x, y, index, "C")

Bond(a::T, b::T) where {T<:Atom2D} = Line2D(a, b)
Bond(a::T, b::T) where {T<:AbstractPoint2D} = Line2D(a, b)
getx(l::Line2D) = l._bx/2
gety(l::Line2D) = l._by/2
getl(l::Line2D) = (l._bx^2 + l._by^2)^0.5

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
include("atoms.jl")
include("bonds.jl")
include("polygons.jl")
include("plotting.jl")

end # module
