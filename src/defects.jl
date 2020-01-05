export make_graphene
export link_relatives!
export make_gatom
export make_gbond
export make_gbond!
export make_bondmatrix
export get_gbond
export make_gpolygon
export make_gpolygon!
export link_bond_polygons!
export make_signature
export make_signature!
# export is_stonewales_bond
# export is_5775_bond
# export is_stonewales
export common_members
export common_relatives
export find_defect
export find_stonewales
export find_5775
export find_divacancy
export find_flower
export find_butterfly
# export make_defect
export filter_relatives_by_type
export filter_members_by_type
export stepout
# export display_dfb

function make_graphene(atom_xy; image_sampling=1, max_bondlength=10.0, frame=1, dataset="dataset")

    indexed_atoms = make_atoms(atom_xy)
    gatom_count = length(indexed_atoms)

    gatoms = make_gatom.(collect(indexed_atoms))

    atom_groups = collect_atom_groups(atom_xy; max_bondlength=max_bondlength)
    bonds, paths = collect_bonds_paths(atom_groups, indexed_atoms)
    if length(bonds) == 0
        graphene = gatoms
        return graphene
    end
    bondmatrix = make_bondmatrix(bonds)
    gbond_count = maximum(bondmatrix)

    gbonds = map(x -> make_gbond!(x, bondmatrix[get_id.(x)...] + gatom_count), [gatoms[collect(y)] for y in Set(sort.(collect.(bonds)))])
    sort!(gbonds, by=get_id)

    polygons = make_polygons(paths)
    polygon_atoms = first.(polygons)
    polygon_bonds = last.(polygons)
    gpolygon_count = length(polygons)
    polygon_gatoms = [gatoms[collect(x)] for x in polygon_atoms]
    polygon_gbonds = [gbonds[x] for x in map(x -> get_gbond(bondmatrix, x), polygon_bonds)]
    polygon_ids = collect(range(1 + gatom_count + gbond_count, length=gpolygon_count))

    gpolygons = map(make_gpolygon!, polygon_gatoms, polygon_gbonds, polygon_ids)

    gpolygon_dict = Dict([Pair(x._id, x) for x in gpolygons])
    map(x -> link_bond_polygons!(gpolygon_dict, x), gbonds)

    graphene = vcat(gatoms, gbonds, gpolygons)

    noa_dict = Dict([Pair(x._id, get_noa(x)) for x in graphene])
    map(x -> make_signature!(x, noa_dict), graphene)

    return graphene
end

function link_relatives!(a, b)
    push!(a._relatives, b._id)
    push!(b._relatives, a._id)
end

function link_relatives!(a, b_vector::Vector)
    link_relatives!.(tuple(a), b_vector)
end

function make_gatom(indexed_atom)
    id = first(indexed_atom)
    atom = last(indexed_atom)
    GAtom(id, getx(atom), gety(atom))
end

function make_gbond(bond_gatoms, gbond_id)
    GBond(gbond_id, mean(get_x.(bond_gatoms)), mean(get_y.(bond_gatoms)), Set(get_id.(bond_gatoms)))
end

function make_gbond!(bond_gatoms::Vector{GAtom}, gbond_id)
    gbond = GBond(gbond_id, mean(get_x.(bond_gatoms)), mean(get_y.(bond_gatoms)))
    link_relatives!(gbond, bond_gatoms)
    return gbond
end

function make_bondmatrix(bonds)
    bond_dict = Dict(map(x -> Pair(x, TupleTools.sort(x)), collect(bonds)))
    bonds_directionless = unique(collect(values(bond_dict)))
    index1 = vcat(first.(bonds_directionless), last.(bonds_directionless))
    index2 = vcat(last.(bonds_directionless), first.(bonds_directionless))
    bondids = collect(1:length(bonds_directionless))
    bondmatrix = sparse(index1, index2, repeat(bondids,2))
end

function get_gbond(bondmatrix, bond_index::Tuple{Int, Int})
    bondmatrix[bond_index...]
end

function get_gbond(bondmatrix, bond_index_vector::Vector)
    [get_gbond(bondmatrix, x) for x in bond_index_vector]
end

function make_gpolygon(polygon_gatoms::Vector{GAtom}, polygon_gbonds::Vector{GBond}, gpolygon_id)
    gpolygon = GPolygon(gpolygon_id, mean(get_x.(polygon_gatoms)), mean(get_y.(polygon_gatoms)), length(polygon_gatoms))
    return gpolygon
end

function make_gpolygon!(polygon_gatoms::Vector{GAtom}, polygon_gbonds::Vector{GBond}, gpolygon_id)
    gpolygon = GPolygon(gpolygon_id, mean(get_x.(polygon_gatoms)), mean(get_y.(polygon_gatoms)), length(polygon_gatoms))
    link_relatives!(gpolygon, polygon_gatoms)
    link_relatives!(gpolygon, polygon_gbonds)
    return gpolygon
end

function link_bond_polygons!(gpolygon_dict, gbond::GBond)
    gpolygon_neighbours = Tuple([gpolygon_dict[x] for x in gbond._relatives if x > gbond._id])
    if length(gpolygon_neighbours) == 2
        link_relatives!(gpolygon_neighbours...)
    end
end

function make_signature(g, noa_dict)
    if isempty(g._relatives)
        signature = "lone atom"
    else
        relatives = [noa_dict[x] for x in g._relatives]
        count_vector = counts(relatives, maximum(relatives))
        signature = join([string(i, "-", count_vector[i], "|") for i in 3:length(count_vector) if !iszero(count_vector[i])])
    end
    return signature
end

function make_signature!(g, noa_dict)
    g._signature = make_signature(g, noa_dict)
end

function find_stonewales(graphene, area)
    defects = GDefect[]
    for g in area
        if get_type(g) == "Bond" && get_signature(g) == "7-2|"
            d = stepout(graphene, g, ["Atom", "Polygon"])
            if sort(get_noa.(filter(ispolygon, graphene[collect(get_members(d))]))) == [5, 5, 7, 7]
                if Set(get_signature.(filter(ispolygon, graphene[collect(get_members(d))]))) == Set(["5-2|6-4|7-1|", "6-3|7-2|"])
                    push!(defects, GDefect(0, get_x(d), get_y(d), get_relatives(d), "", get_frame(d), get_dataset(d), get_noa(d), get_members(d), "Stone-Wales"))
                end
            end
        end
    end
    return defects
end
find_stonewales(graphene) = find_stonewales(graphene, graphene)
# function is_stonewales_bond(graphene, g)
#      condition1 = get_signature(g) == "7-2|"
#      condition2 = get_type(g) == "Bond"
#      condition3 = unique(get_signature.(filter_relatives_by_type(graphene, g, "Atom"))) == ["5-1|7-2|"]
#      condition1 && condition2 && condition3
# end

function find_5775(graphene, area)
    defects = GDefect[]
    for g in area
        if get_type(g) == "Bond" && get_signature(g) == "7-2|"
            gpolygon_members = filter_relatives_by_type(graphene, g, "Polygon")
            if unique(get_signature.(gpolygon_members)) == ["5-1|6-5|7-1|"]
                key_hexagons = filter(x -> get_noa(x) == 6, graphene[collect(get_members(stepout(graphene, g, ["Atom", "Polygon"])))])
                if length(key_hexagons) == 2 && unique(get_signature.(key_hexagons)) == ["5-1|6-3|7-2|"]
                    key_bonds = filter(x -> get_signature(x) == "5-1|7-1|", filter_relatives_by_type(graphene, gpolygon_members, "Bond"))
                    d = stepout(graphene, key_bonds, "Polygon")
                    push!(defects, GDefect(0, get_x(d), get_y(d), get_relatives(d), "", get_frame(d), get_dataset(d), get_noa(d), get_members(d), "5775"))
                end
            end
        end
    end
    return defects
end
find_5775(graphene) = find_5775(graphene, graphene)
# function is_5775_bond(graphene, g)
#      condition1 = get_signature(g) == "7-2|"
#      condition2 = get_type(g) == "Bond"
#      key_gatoms = filter_relatives_by_type(graphene, g, "Atom")
#      condition3 = unique(get_signature.(key_gatoms)) == ["6-1|7-2|"]
#      gpolygon_members = filter_relatives_by_type(graphene, g, "Polygon")
#      condition4 = unique(get_signature.(gpolygon_members)) == ["5-1|6-5|7-1|"]
#      key_gpolygons = filter(x -> get_noa(x) == 6, vcat([filter_relatives_by_type(graphene, x, "Polygon") for x in key_gatoms]...))
#      condition5 = unique(get_signature.(key_gpolygons)) == ["5-1|6-3|7-2|"]
#      condition1 && condition2 && condition3 && condition4 && condition5
# end

function find_divacancy(graphene, area)
    defects = GDefect[]
    for g in area
        if get_type(g) == "Polygon" && get_signature(g) == "5-2|6-6|"
            hexagon_type1 = filter(x -> get_signature(x) == "6-5|8-1|", filter_relatives_by_type(graphene, g, "Polygon"))
            hexagon_type2 = filter(x -> get_signature(x) == "5-1|6-4|8-1|", filter_relatives_by_type(graphene, g, "Polygon"))
            if length(hexagon_type1) == 2 && length(hexagon_type2) == 4
                atoms_type2 = vcat([filter_relatives_by_type(graphene, x, "Atom") for x in hexagon_type2]...)
                if length(atoms_type2) == length(unique(atoms_type2))
                    key_bonds = filter(x -> get_signature(x) == "5-1|8-1|", filter_relatives_by_type(graphene, g, "Bond"))
                    d = stepout(graphene, key_bonds, "Polygon")
                    push!(defects, GDefect(0, get_x(d), get_y(d), get_relatives(d), "", get_frame(d), get_dataset(d), get_noa(d), get_members(d), "Divacancy"))
                end
            end
        end
    end
    return defects
end
find_divacancy(graphene) = find_divacancy(graphene, graphene)

function find_flower(graphene, area)
    defects = GDefect[]
    for g in area
        if get_type(g) == "Atom" && get_signature(g) == "7-3|"
            gpolygon_members = filter_relatives_by_type(graphene, g, "Polygon")
            if unique(get_noa.(gpolygon_members)) == [7] && unique(get_signature.(gpolygon_members)) == ["5-2|6-3|7-2|"]
                d = stepout(graphene, g, ["Bond", "Atom", "Polygon"])
                push!(defects, GDefect(0, get_x(d), get_y(d), get_relatives(d), "", get_frame(d), get_dataset(d), get_noa(d), get_members(d), "Flower"))
            end
        end
    end
    return defects
end
find_flower(graphene) = find_flower(graphene, graphene)

function common_members(graphene, g1, g2, types...)
    m = graphene[collect(intersect(get_members(g1), get_members(g2)))]
    if !isempty(types)
        m = filter(x -> istype(x, types), m)
    end
    return m
end

function common_relatives(graphene, g1, g2, types...)
    r = graphene[collect(intersect(get_relatives(g1), get_relatives(g2)))]
    if !isempty([types...])
        r = filter(x -> istype(x, types...), r)
    end
    return r
end

function find_butterfly(graphene, area)
    defects = GDefect[]
    for g in area
        if get_type(g) == "Polygon" && get_signature(g) == "5-2|7-4|"
            gpolygon_members = filter_relatives_by_type(graphene, g, "Polygon")
            key_polygons = filter(x -> get_noa(x) == 7 && get_signature(x) == "5-2|6-4|7-1|", gpolygon_members)
            if length(key_polygons) == 4
                temp = vcat([common_relatives(graphene, i..., "Atom") for i in subsets(key_polygons, Val(2))]...)
                key_atoms = unique([filter_relatives_by_type(graphene, g, "Atom"); temp])
                d = stepout(graphene, key_atoms, "Polygon")
                push!(defects, GDefect(0, get_x(d), get_y(d), get_relatives(d), "", get_frame(d), get_dataset(d), get_noa(d), get_members(d), "Butterfly"))
            end
        end
    end
    return defects
end
find_butterfly(graphene) = find_butterfly(graphene, graphene)

function find_defect(graphene, type::String)
    id_offset = maximum(get_id.(graphene))
    defects = GDefect[]
    if type == "Stone-Wales"
        defects = [defects; find_stonewales(graphene)]
    elseif type == "5775"
        defects = [defects; find_5775(graphene)]
    elseif type == "Divacancy" || type == "V2(585)"
        defects = [defects; find_divacancy(graphene)]
    elseif type == "Flower" || type == "V2(555-777)"
        defects = [defects; find_flower(graphene)]
    elseif type == "Butterfly" || type == "V2(5555-6-7777)"
        defects = [defects; find_butterfly(graphene)]
    end
    return defects
end
# function find_defect(graphene, type::String)
#     id_offset = maximum(get_id.(graphene))
#     if type == "Stone-Wales"
#         markers = filter(x -> is_stonewales_bond(graphene, x), graphene)
#     elseif type == "5775"
#         markers = filter(x -> is_5775_bond(graphene, x), graphene)
#     else
#         if type == "Flower"
#             type_signature = "7-3|"
#             gtype = "Atom"
#         elseif type == "Butterfly"
#             type_signature = "5-2|7-4|"
#             gtype = "Polygon"
#         elseif type == "Divacancy"
#             type_signature = "5-2|6-6|"
#             gtype = "Polygon"
#         else
#             type_signature = ""
#             gtype = ""
#         end
#         markers = filter(x -> get_signature(x) == type_signature && get_type(x) == gtype, graphene)
#     end
#     map(x -> make_defect(graphene, markers[x], id_offset+x), collect(1:length(markers)))
# end

function find_defect(graphene, types...)
    id_offset = maximum(get_id.(graphene))
    defects = vcat(map(t -> find_defect(graphene, t), types)...)
    for i in 1:length(defects)
        defects[i]._id = id_offset+i
    end
    return defects
end

# function make_defect(graphene, marker, id)
#     if get_type(marker) == "Atom"
#         gbond_members = filter_relatives_by_type(graphene, marker, "Bond")
#         gatom_members = vcat([filter_relatives_by_type(graphene, x, "Atom") for x in gbond_members]...)
#         gpolygon_members = vcat([filter_relatives_by_type(graphene, x, "Polygon") for x in gatom_members]...)
#         if get_signature(marker) == "7-3|"
#             type = "Flower"
#         else
#             type = "Complex"
#         end
#     elseif get_type(marker) == "Bond"
#             gatom_members = filter_relatives_by_type(graphene, marker, "Atom")
#             if get_signature(marker) == "7-2|"
#                 if unique(get_signature.(gatom_members)) == ["5-1|7-2|"]
#                     gpolygon_members = filter_relatives_by_type(graphene, marker, "Polygon")
#                     gpolygon_members = vcat([filter_relatives_by_type(graphene, x, "Polygon") for x in gatom_members]...)
#                     relatives = setdiff(vcat([filter_relatives_by_type(graphene, x, "Polygon") for x in gpolygon_members]...), gpolygon_members)
#                 # if unique(get_noa.(relatives)) == [6]
#                     type = "Stone-Wales"
#                 else
#                     type = "Complex"
#                 end
#                 if unique(get_signature.(gatom_members)) == ["6-1|7-2|"]
#                     gpolygon_members = filter_relatives_by_type(graphene, marker, "Polygon")
#                     gpolygon_members = vcat([filter_relatives_by_type(graphene, x, "Polygon") for x in gpolygon_members]...) |> unique
#                     filter!(x -> get_noa(x) != 6 , gpolygon_members)
#                     relatives = setdiff(vcat([filter_relatives_by_type(graphene, x, "Polygon") for x in gpolygon_members]...), gpolygon_members)
#                     type = "5775"
#                 end
#             else
#                 type = "Complex"
#             end
#     elseif get_type(marker) == "Polygon"
#         if get_signature(marker) == "5-2|7-4|"
#             gatom_members = filter_relatives_by_type(graphene, marker, "Atom")
#             gbond_members = vcat([filter_relatives_by_type(graphene, x, "Bond") for x in gatom_members]...) |> unique
#             gbond_members = vcat(filter(x -> get_signature(x) == "7-2|" , gbond_members), filter_relatives_by_type(graphene, marker, "Bond"))
#             gatom_members = vcat([filter_relatives_by_type(graphene, x, "Atom") for x in gbond_members]...) |> unique
#             gpolygon_members = vcat([filter_relatives_by_type(graphene, x, "Polygon") for x in gatom_members]...) |> unique
#             type = "Butterfly"
#         elseif get_signature(marker) == "5-2|6-6|"
#             gatom_members = filter!(x -> get_signature(x) == "6-2|8-1|", filter_relatives_by_type(graphene, marker, "Atom"))
#             gpolygon_members = vcat(filter_relatives_by_type(graphene, marker, "Polygon"), marker)
#             filter!(x -> get_noa(x) != 6 , gpolygon_members)
#             if length(gatom_members) == 4
#                 points = map(a -> Point2D(get_x(a), get_y(a)), gatom_members)
#                 line_points = map(a -> Point2D(get_x(a), get_y(a)), filter(x -> get_noa(x) == 5, gpolygon_members))
#                 line = Line(line_points...)
#                 points_orientation = map(x -> orientation(line, x), points)
#                 if +(points_orientation...) == 0 || *(points_orientation...) == 1
#                     type = "Divacancy"
#                 else
#                     type = "Complex"
#                 end
#             else
#                 type = "Complex"
#             end
#         else
#             type = "Complex"
#         end
#     else
#         gpolygon_members = filter_relatives_by_type(graphene, marker, "Polygon")
#         type = "Complex"
#     end
#
#     gatom_gbond_members = vcat([filter_relatives_by_type(graphene, x, "Atom", "Bond") for x in gpolygon_members]...) |> unique
#     members = Set(get_id.(vcat(gpolygon_members, gatom_gbond_members)))
#     relatives = setdiff(vcat([filter_relatives_by_type(graphene, x, "Polygon") for x in gpolygon_members]...), gpolygon_members)
#     noa = length(unique(filter(isatom, graphene[collect(members)])))
#     signature = ""
#     if unique(get_noa.(relatives)) != [6]
#         type = "Complex"
#     end
#
#     GDefect(id, get_x(marker), get_y(marker), Set(get_id.(relatives)), signature, get_frame(marker), get_dataset(marker), noa, members, type)
# end

function filter_relatives_by_type(graphene, g, type)
    relatives = graphene[collect(get_relatives(g))]
    filter(x -> get_type(x) == type, relatives)
end
filter_relatives_by_type(graphene, g, types...) = vcat(map(t -> filter_relatives_by_type(graphene, g, t), types)...)
filter_relatives_by_type(graphene, g_vector::Vector, types...) = unique(vcat(map(g -> filter_relatives_by_type(graphene, g, types...), g_vector)...))

function filter_members_by_type(graphene, g, type)
    members = graphene[collect(get_members(g))]
    filter(x -> get_type(x) == type, members)
end
filter_members_by_type(graphene, g, types...) = vcat(map(t -> filter_members_by_type(graphene, g, t), types)...)
filter_members_by_type(graphene, g_vector::Vector, types...) = unique(vcat(map(g -> filter_members_by_type(graphene, g, types...), g_vector)...))

function stepout(graphene, g, n::Int, steps=String[]; id=0)
    gmembers = AbstractGPrimitive[]
    if typeof(g) <: Vector
        members = union(get_members.(g)...)
        relatives = union(get_relatives.(g)...)
        frame = get_frame(first(g))
        dataset = get_dataset(first(g))
    else
        members = get_members(g)
        relatives = get_relatives(g)
        frame = get_frame(g)
        dataset = get_dataset(g)
    end
    if isempty(steps)
        itr = Stateful(cycle(["Atom", "Bond", "Polygon"]))
        if typeof(g) <: AbstractGPrimitive
            for x in itr; x != get_type(g) || break; end
        end
        steps = collect(take(itr, n))
    end

    seed = g
    for t in steps
        seed = filter_relatives_by_type(graphene, seed, t)
        for s in seed
            push!(gmembers, s)
            relatives = union(relatives, get_relatives(s))
        end
    end
    x = mean(get_x.(gmembers))
    y = mean(get_y.(gmembers))
    signature = ""
    if steps[end] != "Atom"
        gmembers = [gmembers; filter(x -> isatom(x)||isbond(x), graphene[collect(relatives)])]
    end
    members = union([get_members.(gmembers); members]...)
    gmembers = [gmembers; filter(x -> isatom(x)||isbond(x), graphene[collect(members)])]
    noa = length(unique(filter(isatom, gmembers)))
    members = Set(get_id.(gmembers))
    relatives = setdiff!(relatives, members)

    g_new = GEntry(id, x, y, relatives, signature, frame, dataset, noa, members)
end
stepout(graphene, g, steps::Vector{String}; id=0) = stepout(graphene, g, 0, steps; id=id)
stepout(graphene, g, steps::String; id=0) = stepout(graphene, g, 0, [steps]; id=id)
stepout(graphene, g; id=0) = stepout(graphene, g, 1; id=id)


# function find_dfb_save(graphene_stack, output_path)
#     divacancy = merge_stack(map(x -> find_isolated_defect(x, "divacancy"), graphene_stack))
#     flower    = merge_stack(map(x -> find_isolated_defect(x, "flower"), graphene_stack))
#     butterfly = merge_stack(map(x -> find_isolated_defect(x, "butterfly"), graphene_stack))
#     dfb       = merge_stack([divacancy, flower, butterfly])
#     FileIO.save(output_path, dfb)
# end
#
# function display_dfb(d_divacancy_merge, d_flower_merge, d_butterfly_merge; shape="")
#     red_line   = zeros(10000)
#     green_line = zeros(10000)
#     blue_line  = zeros(10000)
#     red_countmap   = if ~isempty(d_divacancy_merge) countmap(DataFrame(d_divacancy_merge)[:, :frame]) else [] end
#     green_countmap = if ~isempty(d_flower_merge) countmap(DataFrame(d_flower_merge)[      :, :frame]) else [] end
#     blue_cpoPolygon  = if ~isempty(d_butterfly_merge) countmap(DataFrame(d_butterfly_merge)[:, :frame]) else [] end
#
#     for (k, v) in red_countmap
#             red_line[k] = v
#     end
#     for (k, v) in green_countmap
#             green_line[k] = v
#     end
#     for (k, v) in blue_countmap
#             blue_line[k] = v
#     end
#
#     if shape == "line"
#         red   = red_line
#         green = green_line
#         blue  = blue_line
#     else
#         red   = reshape(red_line, (100, 100))
#         green = reshape(green_line, (100, 100))
#         blue  = reshape(blue_line, (100, 100))
#     end
#
#     rgb_dfb = cat(red, green, blue; dims=3)
#     colorview(RGB, permutedims(rgb_dfb, (3,2,1)))
# end
