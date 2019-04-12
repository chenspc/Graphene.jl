# using Plots
# using Gadfly
# gr()


export plot_atoms, plot_atoms!, plot_bonds, plot_polygons, plot_graphene, make_shapes

# function plot_atoms(graphene_model)
function plot_atoms(atom_xy)
    scatter(atom_xy[1,:], atom_xy[2,:], w=3, xlims=(0,256), ylims=(0,256), aspect_ratio=:equal, leg=false)
end

function plot_atoms!(atom_xy)
    scatter!(atom_xy[1,:], atom_xy[2,:], w=3, xlims=(0,256), ylims=(0,256), aspect_ratio=:equal, leg=false)
end

function plot_bonds(graphene_model)
    return
end

function make_shapes(patoms_collection,indexed_atoms_collection)
    verts_collection = map(x -> index2xy(x, indexed_atoms_collection), patoms_collection)
    shape_collection = map(x -> Shape(x), verts_collection)
end

# function plot_polygons(graphene_model)
function plot_polygons(patoms_collection, indexed_atoms_collection)
    shape_collection = make_shapes(patoms_collection,indexed_atoms_collection)
    plot([],[], xlim =(0,256), ylim=(0,256), aspect_ratio=:equal,leg=false)
    plot!(shape_collection, opacity=.5)
end

function plot_graphene(atom_xy, patoms_collection, indexed_atoms_collection)
    # atom_xy = map(xx -> (getx(xx), gety(xx)), values(indexed_atoms_collection))
    plot_polygons(patoms_collection, indexed_atoms_collection)
    plot_atoms!(atom_xy)
end
