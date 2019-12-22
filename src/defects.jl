# export make_graphene
export index2xy
# export pack_relatives
# export unpack_relatives
# export generate_signature
# export find_defect
# export merge_stack
# export find_isolated_defect
# export isolated_flower
# export isolated_butterfly
# export display_dfb

function make_graphene(atom_xy; image_sampling=1, max_bondlength=10.0, frame=1, dataset="standalone_dataset")

    indexed_atoms = make_atoms(atom_xy)

    atom_groups = collect_atom_groups(atom_xy; max_bondlength=max_bondlength)
    bonds, paths = collect_bonds_paths(atom_groups, indexed_atoms)
    polygons = make_polygons(paths)

    polygon_atoms = first.(polygons)
    polygon_bonds = last.(polygons)

    polygon_atomsxy = map(x->index2xy(x, indexed_atoms), polygon_atoms)
    polygon_x       = map(mean, map(x -> first.(x), polygon_atomsxy))
    polygon_y       = map(mean, map(x -> last.(x), polygon_atomsxy))

    bonds_directional = collect(Iterators.flatten(polygon_bonds))
    f = x -> x=>Tuple(sort(collect(x)))
    bond_dict = Dict(map(f, bonds_directional))
    all_bonds = unique(collect(values(bond_dict)))

    numberof_atomids    = length(indexed_atoms)
    numberof_bondids    = length(all_bonds)
    numberof_polygonids = length(polygons)

    atomid_dict   = Dict(map(x -> UInt32(x) => x, 1:numberof_atomids))
    atomtype_dict = Dict(map(x -> UInt32(x) => "Atom", 1:numberof_atomids))
    atomx_dict    = Dict(map(x -> UInt32(x) => getx(indexed_atoms[x]), 1:numberof_atomids))
    atomy_dict    = Dict(map(x -> UInt32(x) => gety(indexed_atoms[x]), 1:numberof_atomids))

    bondid_dict   = Dict(map(x -> UInt32(numberof_atomids + x) => all_bonds[x], 1:numberof_bondids))
    bondtype_dict = Dict(map(x -> UInt32(numberof_atomids + x) => "Bond", 1:numberof_bondids))
    bondx_dict    = Dict(map(x -> UInt32(numberof_atomids + x) => getx(Bond(indexed_atoms[first(all_bonds[x])], indexed_atoms[last(all_bonds[x])])), 1:numberof_bondids))
    bondy_dict    = Dict(map(x -> UInt32(numberof_atomids + x) => gety(Bond(indexed_atoms[first(all_bonds[x])], indexed_atoms[last(all_bonds[x])])), 1:numberof_bondids))

    polygonid_dict   = Dict(map(x -> UInt32(numberof_atomids + numberof_bondids + x) => polygons[x], 1:numberof_polygonids))
    polygontype_dict = Dict(map(x -> UInt32(numberof_atomids + numberof_bondids + x) => "Polygon", 1:numberof_polygonids))
    polygonx_dict    = Dict(map(x -> UInt32(numberof_atomids + numberof_bondids + x) => polygon_x[x], 1:numberof_polygonids))
    polygony_dict    = Dict(map(x -> UInt32(numberof_atomids + numberof_bondids + x) => polygon_y[x], 1:numberof_polygonids))

    atomid_dict_reversed = Dict( v => k for (k, v) in atomid_dict)
    bondid_dict_reversed = merge(Dict(v => k for (k, v) in bondid_dict),
                                 Dict(reverse(v) => k for (k, v) in bondid_dict))

    a2b_dict = Dict{UInt32, Vector{UInt32}}()
    for (k, v) in bondid_dict
        for vi in v
            if haskey(a2b_dict, atomid_dict_reversed[vi])
                push!(a2b_dict[atomid_dict_reversed[vi]],k)
            else
                a2b_dict[atomid_dict_reversed[vi]] = [k]
            end
        end
    end

    a2p_dict = Dict{UInt32, Vector{UInt32}}()
    for (k, v) in polygonid_dict
        for vi in first(v)
            if haskey(a2p_dict, atomid_dict_reversed[vi])
                # push!(a2p_dict[vi],k)
                push!(a2p_dict[atomid_dict_reversed[vi]],k)
            else
                # a2p_dict[vi] = [k]
                a2p_dict[atomid_dict_reversed[vi]] = [k]
            end
        end
    end

    a2bp_dict  = Dict{UInt32, Vector{UInt32}}()
    a2noa_dict = Dict{UInt32, Int}()
    for (k, v) in a2b_dict
        if haskey(a2p_dict, k)
            a2bp_dict[k] = vcat(a2b_dict[k], a2p_dict[k])
        else
            a2bp_dict[k] = a2b_dict[k]
        end
        a2noa_dict[k] = 1
    end

    b2p_dict = Dict{UInt32, Vector{UInt32}}()
    for (k, v) in polygonid_dict
        for vi in last(v)
            if haskey(b2p_dict, bondid_dict_reversed[bond_dict[vi]])
                push!(b2p_dict[bondid_dict_reversed[bond_dict[vi]]],k)
            else
                b2p_dict[bondid_dict_reversed[bond_dict[vi]]] = [k]
            end
        end
    end

    b2a_dict = Dict{UInt32, Vector{UInt32}}()
    for (k, v) in bondid_dict
        b2a_dict[k] = [atomid_dict_reversed[first(v)], atomid_dict_reversed[last(v)]]
    end

    b2ap_dict  = Dict{UInt32, Vector{UInt32}}()
    b2noa_dict = Dict{UInt32, Int}()
    for (k, v) in b2a_dict
        b2ap_dict[k] = vcat(b2a_dict[k], b2p_dict[k])
        b2noa_dict[k] = 2
    end


    p2ab_dict  = Dict{UInt32, Vector{UInt32}}()
    p2noa_dict = Dict{UInt32, Int}()
    for (k, v) in polygonid_dict
        p2ab_dict[k] = vcat(map(x -> atomid_dict_reversed[x], first(v)), map(x -> bondid_dict_reversed[bond_dict[x]], last(v)))
        p2noa_dict[k] = length(first(v))
    end

    p2p_dict = Dict{UInt32, Vector{UInt32}}()
    for (k, v) in p2ab_dict
        p2p_dict[k] = unique!(filter!(x -> x != 0 && x != k,collect(Iterators.flatten(map(x -> get(b2p_dict, x, 0), v)))))
    end

    p2abp_dict = Dict{UInt32, Vector{UInt32}}()
    for (k, v) in p2ab_dict
        p2abp_dict[k] = vcat(p2ab_dict[k], p2p_dict[k])
    end

    g2type_dict = merge(atomtype_dict, bondtype_dict, polygontype_dict)
    g2x_dict    = merge(atomx_dict, bondx_dict, polygonx_dict)
    g2y_dict    = merge(atomy_dict, bondy_dict, polygony_dict)
    idkeys = collect(keys(g2x_dict))
    g2relatives_dict_temp = merge(a2bp_dict, b2ap_dict, p2abp_dict)
    g2relatives_dict = Dict{UInt32, String}()
    for (k, v) in g2relatives_dict_temp
        g2relatives_dict[k] = pack_relatives(v)
    end
    g2noa_dict = merge(a2noa_dict, b2noa_dict, p2noa_dict)

    g2p_dict = merge(a2p_dict, b2p_dict, p2p_dict)
    g2signature_dict = Dict{UInt32, String}()
    max_pside  = 21
    max_pcount = 255
    empty_signature = fill(UInt8(0), max_pside)

    for k in idkeys
        if haskey(g2p_dict, k)
            pcount = countmap(values(map(x -> g2noa_dict[x], g2p_dict[k])))
            g2signature_dict[k] = generate_signature(pcount; max_pside=max_pside, max_pcount=max_pcount)
        end
    end

    t = table((dataset=fill(dataset, length(idkeys)), frame=fill(frame, length(idkeys)), id=idkeys, type=sort_dict(g2type_dict, idkeys), noa=sort_dict(g2noa_dict, idkeys), x=sort_dict(g2x_dict, idkeys), y=sort_dict(g2y_dict, idkeys), relatives=sort_dict(g2relatives_dict, idkeys), signature=sort_dict(g2signature_dict, idkeys)); pkey = :id)
    return t
end


function index2xy(x, indexed_atoms)
    if length(x) == 1
        atom   = indexed_atoms[x[1]]
        output = [(getx(atom), gety(atom))]
    else
        atoms  = indexed_atoms[x]
        output = map(xx -> (getx(xx), gety(xx)), atoms)
    end
    return output
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
