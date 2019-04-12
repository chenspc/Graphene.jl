using Graphene
using Test
using BenchmarkTools
# using Gadfly
using Plots
# using Luxor
# plotly()

# map(getx, values(indexed_atoms_collection))
# map(gety, values(indexed_atoms_collection))

# scatter(atom_xy[1,:], atom_xy[2,:], w=3, xlims=(0,256), ylims=(0,256), aspect_ratio=:equal)

patoms_collection
pbonds_collection

plot_atoms(atom_xy)
polygon_collection[1][1]
polygon = patoms_collection[1]
verts = index2xy(polygon, indexed_atoms_collection)

# vertsx = map(first, verts)
# vertsy = map(last, verts)
# push!(vertsx, first(vertsx), NaN)
# @benchmark Gadfly.plot(x=vertsx, y=vertsy, Coord.cartesian(xmin=0, xmax=256, ymin=0, ymax=256, fixed=true),Geom.polygon)

polygon = patoms_collection[1]
indexed_atoms_collection[polygon][1]
length([indexed_atoms_collection[polygon][1]])

# verts = map(x -> index2xy(x, indexed_atoms_collection), polygon)
verts = []
verts_new = index2xy(polygon, indexed_atoms_collection)
@test_broken verts == verts_new

# @benchmark plot_polygons(patoms_collection, indexed_atoms_collection)

plot_graphene(atom_xy, patoms_collection, indexed_atoms_collection)
