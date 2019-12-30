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
export find_defect
export make_defect
export filter_relatives_by_type
# export merge_stack
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
        signature = join([string(i, "-", count_vector[i], "|") for i in 1:length(count_vector) if !iszero(count_vector[i])])
    end
    return signature
end

function make_signature!(g, noa_dict)
    g._signature = make_signature(g, noa_dict)
end

function find_defect(graphene, type::String)
    id_offset = maximum(get_id.(graphene))
    if type == "Flower"
        type_signature = "2-3|7-3|"
        gtype = "Atom"
    elseif type == "Butterfly"
        type_signature = "1-6|2-6|5-2|7-4|"
        gtype = "Polygon"
    elseif type == "Divacancy"
        type_signature = "1-8|2-8|5-2|6-6|"
        gtype = "Polygon"
    # elseif type == "S-W"
    #     type_signature = "2-3|5-1|7-2|"
    #     gtype = "Bond"
    # elseif type == "5775"
    #     type_signature = "2-3|6-1|7-2|"
    #     gtype = "Bond"
    else
        type_signature = ""
        gtype = ""
    end
    markers = filter(x -> get_signature(x) == type_signature && get_type(x) == gtype, graphene)
    map(x -> make_defect(graphene, markers[x], id_offset+x), collect(1:length(markers)))
end

function make_defect(graphene, marker, id)
    if get_type(marker) == "Atom"
        gbond_members = filter_relatives_by_type(graphene, marker, "Bond")
        gatom_members = vcat([filter_relatives_by_type(graphene, x, "Atom") for x in gbond_members]...)
        gpolygon_members = vcat([filter_relatives_by_type(graphene, x, "Polygon") for x in gatom_members]...)
        if get_signature(marker) == "2-3|7-3|"
            type = "Flower"
        else
            type = "Complex"
        end
    # elseif get_type(marker) == "Bond"
    #     gatom_members = filter_relatives_by_type(graphene, marker, "Atom")
    #     gpolygon_members = vcat([filter_relatives_by_type(graphene, x, "Polygon") for x in gatom_members]...)
    #     if get_signature(marker) == "1-2|7-2|"
    #         if unique(get_signature.(gatom_members)) == "2-3|5-1|7-2|"
    #             type = "S-W"
    #         elseif unique(get_signature.(gatom_members)) == "2-3|6-1|7-2|"
    #             type = "5775"
    #         else
    #             type = "Complex"
    #         end
    #     end
    elseif get_type(marker) == "Polygon"
        if get_signature(marker) == "1-6|2-6|5-2|7-4|"
            gatom_members = filter_relatives_by_type(graphene, marker, "Atom")
            gbond_members = vcat([filter_relatives_by_type(graphene, x, "Bond") for x in gatom_members]...) |> unique
            gbond_members = vcat(filter(x -> get_signature(x) == "1-2|7-2|" , gbond_members), filter_relatives_by_type(graphene, marker, "Bond"))
            gatom_members = vcat([filter_relatives_by_type(graphene, x, "Atom") for x in gbond_members]...) |> unique
            gpolygon_members = vcat([filter_relatives_by_type(graphene, x, "Polygon") for x in gatom_members]...) |> unique
            type = "Butterfly"
        elseif get_signature(marker) == "1-8|2-8|5-2|6-6|"
            gatom_members = filter!(x -> get_signature(x) == "2-3|6-2|8-1|", filter_relatives_by_type(graphene, marker, "Atom"))
            gpolygon_members = vcat(filter_relatives_by_type(graphene, marker, "Polygon"), marker)
            filter!(x -> get_noa(x) != 6 , gpolygon_members)
            if length(gatom_members) == 4
                points = map(a -> Point2D(get_x(a), get_y(a)), gatom_members)
                line_points = map(a -> Point2D(get_x(a), get_y(a)), filter(x -> get_noa(x) == 5, gpolygon_members))
                line = Line(line_points...)
                points_orientation = map(x -> orientation(line, x), points)
                if +(points_orientation...) == 0 || *(points_orientation...) == 1
                    type = "Divacancy"
                else
                    type = "Complex"
                end
            else
                type = "Complex"
            end
        else
            type = "Complex"
        end
    else
        gpolygon_members = filter_relatives_by_type(graphene, marker, "Polygon")
        type = "Complex"
    end

    gatom_gbond_members = vcat([filter_relatives_by_type(graphene, x, "Atom", "Bond") for x in gpolygon_members]...) |> unique
    members = Set(get_id.(vcat(gpolygon_members, gatom_gbond_members)))
    relatives = setdiff(vcat([filter_relatives_by_type(graphene, x, "Polygon") for x in gpolygon_members]...), gpolygon_members)
    noa = length(unique(filter(isatom, graphene[collect(members)])))
    signature = ""
    if unique(get_noa.(relatives)) != [6]
        type = "Complex"
    end

    GDefect(id, get_x(marker), get_y(marker), Set(get_id.(relatives)), signature, get_frame(marker), get_dataset(marker), noa, members, type)
end

function filter_relatives_by_type(graphene, g, type)
    relatives = graphene[collect(get_relatives(g))]
    filter(x -> get_type(x) == type, relatives)
end

function filter_relatives_by_type(graphene, g, types...)
    vcat(map(t -> filter_relatives_by_type(graphene, g, t), types)...)
end

# function merge_stack(t_stack)
#     stack_length = length(t_stack)
#     if stack_length > 1
#         t_merge = t_stack[1]
#         for i in 2:stack_length
#             t_merge = merge(t_merge, t_stack[i])
#         end
#     else
#         t_merge = t_stack
#     end
#     return t_merge
# end
#
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
