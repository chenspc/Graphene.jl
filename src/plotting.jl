export plot_atom
export plot_atom!
export plot_graphene

@recipe function f(g::Vector{T}) where T<:AbstractGPrimitive
         get_x.(g), get_y.(g)
         end

function plot_atom!(plt, g::Vector{T}; kw...) where T<:AbstractGPrimitive
    atoms = filter(isatom, g)
    bonds = filter(isbond, g)
    polygons = filter(ispolygon, g)
    scatter!(plt, get_x.(atoms), get_y.(atoms); kw...)
end
plot_atom!(plt, g::T; kw...) where T<:AbstractGPrimitive = plot_atom!(plt, [g], kw...)

function plot_atom!(g::Vector{T}; kw...) where T<:AbstractGPrimitive
    atoms = filter(isatom, g)
    bonds = filter(isbond, g)
    polygons = filter(ispolygon, g)
    scatter!(get_x.(atoms), get_y.(atoms); kw...)
end
plot_atom!(g::T; kw...) where T<:AbstractGPrimitive = plot_atom!([g], kw...)

function plot_atom(g::Vector{T}; kw...) where T<:AbstractGPrimitive
    atoms = filter(isatom, g)
    bonds = filter(isbond, g)
    polygons = filter(ispolygon, g)
    scatter(get_x.(atoms), get_y.(atoms); kw...)
end
plot_atom(g::T; kw...) where T<:AbstractGPrimitive = plot_atom([g], kw...)

function plot_graphene(graphene, g::Vector{T}; kw...) where T<:AbstractGPrimitive
    plt = plot([], aspect_ratio = 1, xlims=(0,256), ylims=(0,256), framestyle=:box, legend=false)

    g_palette = [
    colorant"#6E6E6E", # 1
    colorant"#DF0101", # 2
    colorant"#FBD606", # 3
    colorant"#BF2626", # 4
    colorant"#FEE081", # 5
    colorant"#B3BED2", # 6
    colorant"#353D52", # 7
    colorant"#0D621E", # 8
    colorant"#934BB2", # 9
    colorant"#3EAE7E", # 10
    colorant"#EFBAA0", # 11
    colorant"#602E56", # 12
    colorant"#3f9778", # 13
    colorant"#FF8000", # 14
    colorant"#5FB404", # 15
    colorant"#182854", # 16
    colorant"#ACAEB5", # 17
    colorant"#5A7460", # 18
    colorant"#B93F01", # 19
    colorant"#1897AB", # 20
    colorant"#FFFFFF"] # >20

    gatoms = filter(isatom, g)
    # gbonds = filter(isbond, g)
    gpolygons = sort(filter(ispolygon, g), rev=true)
    for n in gpolygons
        points = [[get_x(p), get_y(p)] for p in filter_relatives_by_type(graphene, n, :atom)]
        hull = convex_hull(points)
        noa = get_noa(n)
        plot!(plt, VPolygon(hull), alpha=1, color=g_palette[min(noa,21)])
    end
    if !isempty(gatoms)
        plot_atom!(plt, gatoms, color=g_palette[1])
    end
    return plt
end
plot_graphene(graphene, g::T; kw...) where T<:AbstractGEntry = plot_polygon(graphene, [g]; kw...)

function plot_graphene(graphene, g_vector::Vector{GDefect}; kw...)
    temp = AbstractGPrimitive[]
    for g in g_vector
        temp = [graphene[collect(get_members(g))]; temp]
    end
    plot_graphene(graphene, temp; kw...)
end
:atom
