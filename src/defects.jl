export make_graphene
export link_relatives!
export make_gatom
export make_gbond
export make_gbond!
export make_bondmatrix
export get_gbond
export make_gpolygon
export make_gpolygon!
# export pack_relatives
# export unpack_relatives
# export generate_signature
# export find_defect
# export merge_stack
# export find_isolated_defect
# export isolated_flower
# export isolated_butterfly
# export display_dfb

function make_graphene(atom_xy; image_sampling=1, max_bondlength=10.0, frame=1, dataset="dataset")

    indexed_atoms = make_atoms(atom_xy)
    numberof_gatoms = length(indexed_atoms)

    gatoms = make_gatom.(collect(indexed_atoms))

    atom_groups = collect_atom_groups(atom_xy; max_bondlength=max_bondlength)
    bonds, paths = collect_bonds_paths(atom_groups, indexed_atoms)
    bondmatrix = make_bondmatrix(bonds)
    numberof_gbonds = maximum(bondmatrix)

    gbonds = map(x -> make_gbond!(x, bondmatrix[get_id.(x)...] + numberof_gatoms), [gatoms[collect(y)] for y in bonds if first(y) < last(y)])

    polygons = make_polygons(paths)
    polygon_atoms = first.(polygons)
    polygon_bonds = last.(polygons)
    numberof_gpolygons = length(polygons)
    polygon_gatoms = [gatoms[collect(x)] for x in polygon_atoms]
    polygon_gbonds = [gbonds[x] for x in map(x -> get_gbond(bondmatrix, x), polygon_bonds)]
    polygon_ids = collect(range(1+numberof_gatoms+numberof_gbonds, length=numberof_gpolygons))

    gpolygons = map(make_gpolygon!, polygon_gatoms, polygon_gbonds, polygon_ids)

    # g2signature_dict = Dict{UInt32, String}()
    # max_pside  = 21
    # max_pcount = 255
    # empty_signature = fill(UInt8(0), max_pside)
    #
    # for k in idkeys
    #     if haskey(g2p_dict, k)
    #         pcount = countmap(values(map(x -> g2noa_dict[x], g2p_dict[k])))
    #         g2signature_dict[k] = generate_signature(pcount; max_pside=max_pside, max_pcount=max_pcount)
    #     end
    # end
    #
    # t = table((dataset=fill(dataset, length(idkeys)), frame=fill(frame, length(idkeys)), id=idkeys, type=sort_dict(g2type_dict, idkeys), noa=sort_dict(g2noa_dict, idkeys), x=sort_dict(g2x_dict, idkeys), y=sort_dict(g2y_dict, idkeys), relatives=sort_dict(g2relatives_dict, idkeys), signature=sort_dict(g2signature_dict, idkeys)); pkey = :id)
    # return t
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

function make_signature(g)
    body
end

function make_signature!(g)
    signature = make_signature(g)
    g._signature = signature
    return signature
end


# function sort_dict(dict, idkeys)
#     map(x -> get(dict, x, 0), idkeys)
# end
#
# function pack_relatives(relatives::Array{UInt32,1})
#     base64encode(relatives)
# end
#
# function unpack_relatives(relatives_str::String)
#     reinterpret(UInt32, base64decode(relatives_str)) |> Vector
# end
#
# function generate_signature(pcount; max_pside=21, max_pcount=255)
#     signature_str = fill(UInt8(0), max_pside)
#     for (ks, vs) in pcount
#         signature_str[min(ks, max_pside)] = min(vs, max_pcount)
#     end
#     encoded_signature_str = base64encode(signature_str)
#     return encoded_signature_str
# end
#
# function find_defect(g_table::IndexedTable, type::String)
#     if type == "flower"
#         type_signature = generate_signature(Dict(7 => 3))
#         defects = filter(x -> (x.signature == type_signature), g_table)
#     elseif type == "butterfly"
#         type_signature = generate_signature(Dict(7 => 4, 5 => 2))
#         defects = filter(x -> (x.signature == type_signature), g_table)
#     elseif type == "both"
#         type_signature1 = generate_signature(Dict(7 => 3))
#         type_signature2 = generate_signature(Dict(7 => 4, 5 => 2))
#         defects = filter(x -> (x.signature == type_signature1) || (x.signature == type_signature2), g_table)
#     elseif type == "divacancy"
#         type_signature = generate_signature(Dict(6 => 10))
#         defects = filter(x -> (x.signature == type_signature) && (x.noa == 14), g_table)
#     elseif type == "all"
#         type_signature1 = generate_signature(Dict(7 => 3))
#         type_signature2 = generate_signature(Dict(7 => 4, 5 => 2))
#         type_signature3 = generate_signature(Dict(6 => 10))
#         defects = filter(x -> (x.signature == type_signature1) || (x.signature == type_signature2) || ((x.signature == type_signature3) && (x.noa == 14)), g_table)
#     else
#         type_signature = type
#         defects = filter(x -> (x.signature == type_signature), g_table)
#     end
#     return defects
# end
#
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
# function find_isolated_defect(g_table::IndexedTable, type::String)
#     g_df = DataFrame(g_table)
#     if type == "flower"
#         type_signature = generate_signature(Dict(7 => 3))
#         defects = filter(x -> (x.signature == type_signature), g_table)
#         isolated_defects = filter(x -> isolated_flower(g_df, x), defects)
#     elseif type == "butterfly"
#         type_signature = generate_signature(Dict(7 => 4, 5 => 2))
#         defects = filter(x -> (x.signature == type_signature), g_table)
#         isolated_defects = filter(x -> isolated_butterfly(g_df, x), defects)
#     elseif type == "both"
#         # Not yet modified for isolated defects
#         type_signature1 = generate_signature(Dict(7 => 3))
#         type_signature2 = generate_signature(Dict(7 => 4, 5 => 2))
#         defects = filter(x -> (x.signature == type_signature1) || (x.signature == type_signature2), g_table)
#     elseif type == "divacancy"
#         type_signature1 = generate_signature(Dict(6 => 10))
#         type_signature2 = generate_signature(Dict(5 => 2, 6 => 6))
#         defects = filter(x -> (x.signature == type_signature1) && (x.noa == 14) || (x.signature == type_signature2) && (x.noa == 8) , g_table)
#         isolated_defects = filter(x -> isolated_divacancy(g_df, x), defects)
#     elseif type == "all"
#         # Not yet modified for isolated defects
#         type_signature1 = generate_signature(Dict(7 => 3))
#         type_signature2 = generate_signature(Dict(7 => 4, 5 => 2))
#         type_signature3 = generate_signature(Dict(6 => 10))
#         type_signature4 = generate_signature(Dict(5 => 2, 6 => 6))
#         defects = filter(x -> (x.signature == type_signature1) || (x.signature == type_signature2) || ((x.signature == type_signature3) && (x.noa == 14)) || ((x.signature == type_signature4) && (x.noa == 8)), g_table)
#     else
#         type_signature = type
#         defects = filter(x -> (x.signature == type_signature), g_table)
#     end
#     return isolated_defects
# end
#
# function isolated_flower(g_df::DataFrame, defect)
#     is_isolated = false
#
#     relatives = unpack_relatives(defect[:relatives])
#     relatives_df = g_df[findall(in(relatives), g_df.id), :]
#     selected_relatives_df = relatives_df[findall(in(7), relatives_df.noa), :]
#
#     if unique(selected_relatives_df[:signature]) == [generate_signature(Dict(5 => 2, 6 => 3, 7 => 2))]
#         all_relatives = unique(collect(Iterators.flatten(collect.(map(x -> unpack_relatives(x), selected_relatives_df[:relatives])))))
#         all_relatives_df = g_df[findall(in(all_relatives), g_df.id), :]
#         pentagon_relatives_df = all_relatives_df[findall(in(5), all_relatives_df.noa), :]
#         if nrow(pentagon_relatives_df) == 3
#             is_isolated = true
#         end
#     end
#
#     is_isolated
# end
#
# function isolated_butterfly(g_df::DataFrame, defect)
#     is_isolated = false
#
#     relatives = unpack_relatives(defect[:relatives])
#     relatives_df = g_df[findall(in(relatives), g_df.id), :]
#
#     all_relatives = unique(collect(Iterators.flatten(collect.(map(x -> unpack_relatives(x), relatives_df[:relatives])))))
#     all_relatives_df = g_df[findall(in(all_relatives), g_df.id), :]
#
#     pentagon_relatives_df = all_relatives_df[findall(in(5), all_relatives_df.noa), :]
#     heptagon_relatives_df = all_relatives_df[findall(in(7), all_relatives_df.noa), :]
#
#     pentagons_sigature = pentagon_relatives_df[:signature]
#     heptagons_sigature = heptagon_relatives_df[:signature]
#
#     pentagons_check = countmap([c for c in pentagons_sigature]) == Dict(generate_signature(Dict(6 => 3, 7 => 2)) => 4)
#     heptagons_check = countmap([c for c in heptagons_sigature]) == Dict(generate_signature(Dict(5 => 2, 6 => 4, 7 => 1)) => 4)
#
#     if pentagons_check && heptagons_check
#         is_isolated = true
#     end
#
#     is_isolated
# end
#
# function isolated_divacancy(g_df::DataFrame, defect)
#     is_isolated = false
#
#     if defect[:noa] == 14
#         is_isolated = true
#     else
#         relatives = unpack_relatives(defect[:relatives])
#         relatives_df = g_df[findall(in(relatives), g_df.id), :]
#         selected_relatives_df = relatives_df[findall(in(5), relatives_df.noa), :]
#
#         if unique(selected_relatives_df[:signature]) == [generate_signature(Dict( 6 => 4, 8 => 1))]
#             is_isolated = true
#         end
#     end
#
#     is_isolated
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
#     blue_countmap  = if ~isempty(d_butterfly_merge) countmap(DataFrame(d_butterfly_merge)[:, :frame]) else [] end
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
