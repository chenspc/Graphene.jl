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
# using Distances
# using Images
using BenchmarkTools

export CAtom,
       CBond,
       CPolygon

# Change this to enable debugging
const DEBUG = false

include("data_inport.jl")
include("id_generator.jl")
include("atoms.jl")
include("bonds.jl")
include("polygons.jl")

end # module
