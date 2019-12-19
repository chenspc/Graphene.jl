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

    indexed_atoms_collection = make_atoms(atom_xy)

    atom_groups_collection = collect_atom_groups(atom_xy; max_bondlength=max_bondlength)
    bonds_collection, paths_collection = collect_bonds_paths(atom_groups_collection, indexed_atoms_collection)
    polygon_collection = collect_polygons(paths_collection)

    patoms_collection = first.(polygon_collection)
    pbonds_collection = last.(polygon_collection)

    patomsxy_collection = map(x->index2xy(x, indexed_atoms_collection), patoms_collection)
    px_collection = map(mean, map(x -> first.(x), patomsxy_collection))
    py_collection = map(mean, map(x -> last.(x), patomsxy_collection))

    all_bonds_directional = collect(Iterators.flatten(pbonds_collection))
    f = x -> x=>Tuple(sort(collect(x)))
    bond_dict = Dict(map(f, all_bonds_directional))
    all_bonds = unique(collect(values(bond_dict)))

    numberof_atomids = length(indexed_atoms_collection)
    numberof_bondids = length(all_bonds)
    numberof_polygonids = length(polygon_collection)

    atomid_dict = Dict(map(x -> UInt32(x) => x, 1:numberof_atomids))
    atomtype_dict = Dict(map(x -> UInt32(x) => "Atom", 1:numberof_atomids))
    atomx_dict = Dict(map(x -> UInt32(x) => getx(indexed_atoms_collection[x]), 1:numberof_atomids))
    atomy_dict = Dict(map(x -> UInt32(x) => gety(indexed_atoms_collection[x]), 1:numberof_atomids))

    bondid_dict = Dict(map(x -> UInt32(numberof_atomids + x) => all_bonds[x], 1:numberof_bondids))
    bondtype_dict = Dict(map(x -> UInt32(numberof_atomids + x) => "Bond", 1:numberof_bondids))
    bondx_dict = Dict(map(x -> UInt32(numberof_atomids + x) => getx(Bond(indexed_atoms_collection[first(all_bonds[x])], indexed_atoms_collection[last(all_bonds[x])])), 1:numberof_bondids))
    bondy_dict = Dict(map(x -> UInt32(numberof_atomids + x) => gety(Bond(indexed_atoms_collection[first(all_bonds[x])], indexed_atoms_collection[last(all_bonds[x])])), 1:numberof_bondids))

    polygonid_dict = Dict(map(x -> UInt32(numberof_atomids + numberof_bondids + x) => polygon_collection[x], 1:numberof_polygonids))
    polygontype_dict = Dict(map(x -> UInt32(numberof_atomids + numberof_bondids + x) => "Polygon", 1:numberof_polygonids))
    polygonx_dict = Dict(map(x -> UInt32(numberof_atomids + numberof_bondids + x) => px_collection[x], 1:numberof_polygonids))
    polygony_dict = Dict(map(x -> UInt32(numberof_atomids + numberof_bondids + x) => py_collection[x], 1:numberof_polygonids))

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

    a2bp_dict = Dict{UInt32, Vector{UInt32}}()
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

    b2ap_dict = Dict{UInt32, Vector{UInt32}}()
    b2noa_dict = Dict{UInt32, Int}()
    for (k, v) in b2a_dict
        b2ap_dict[k] = vcat(b2a_dict[k], b2p_dict[k])
        b2noa_dict[k] = 2
    end


    p2ab_dict = Dict{UInt32, Vector{UInt32}}()
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
    g2x_dict = merge(atomx_dict, bondx_dict, polygonx_dict)
    g2y_dict = merge(atomy_dict, bondy_dict, polygony_dict)
    idkeys = collect(keys(g2x_dict))
    g2relatives_dict_temp = merge(a2bp_dict, b2ap_dict, p2abp_dict)
    g2relatives_dict = Dict{UInt32, String}()
    for (k, v) in g2relatives_dict_temp
        g2relatives_dict[k] = pack_relatives(v)
    end
    g2noa_dict = merge(a2noa_dict, b2noa_dict, p2noa_dict)

    g2p_dict = merge(a2p_dict, b2p_dict, p2p_dict)
    g2signature_dict = Dict{UInt32, String}()
    max_pside = 21
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


function index2xy(x, indexed_atoms_collection)
    if length(x) == 1
        atom = indexed_atoms_collection[x[1]]
        output = [(getx(atom), gety(atom))]
    else
        atoms = indexed_atoms_collection[x]
        output = map(xx -> (getx(xx), gety(xx)), atoms)
    end
    return output
end
