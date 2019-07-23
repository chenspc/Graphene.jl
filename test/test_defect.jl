using Graphene
using Test
using BenchmarkTools
using Profile
using JuliaDB


idkeys = g._id
g2type_dict = g._type
g2x_dict = g._x
g2y_dict = g._y
g2relatives_dict = g._relatives
g2signature_dict = g._signature
g2noa_dict = g._noa
frame = g._frame
dataset = g._dataset

Set(idkeys) == Set(keys(g2type_dict))
Set(idkeys) == Set(keys(g2x_dict))
Set(idkeys) == Set(keys(g2y_dict))
Set(idkeys) == Set(keys(g2relatives_dict))
Set(idkeys) == Set(keys(g2signature_dict))
Set(idkeys) == Set(keys(g2noa_dict))

for i in idkeys

end

map(, idkeys)

map(x -> get(g2type_dict, x, 0), idkeys)
function sort_dict(dict, idkeys)
    map(x -> get(dict, x, 0), idkeys)
end

sort_dict(g2type_dict, idkeys)

t = table((id=idkeys, type=sort_dict(g2type_dict, idkeys), noa=sort_dict(g2noa_dict, idkeys), x=sort_dict(g2x_dict, idkeys), y=sort_dict(g2y_dict, idkeys), relatives=sort_dict(g2relatives_dict, idkeys), signature=sort_dict(g2signature_dict, idkeys), frame=fill(frame, length(idkeys)), dataset=fill(dataset, length(idkeys))); pkey = :id)
t[1]

make_graphene()
@benchmark t = make_graphene(data)


length(idkeys)
fill(1, size(idkeys, 1))
fill("", size(idkeys, 1))

A = idkeys
B = collect(keys(g2type_dict))
filter(x -> !(x in A), B)

g2type_dict[0x25d1f8e46df6eb03]
g2type_dict[0xbadd3b54ef02de11]

g2x_dict[0x25d1f8e46df6eb03]
g2x_dict[0xbadd3b54ef02de11]
g2y_dict[0x25d1f8e46df6eb03]
g2y_dict[0xbadd3b54ef02de11]
g2x_dict
g2y_dict

g2relatives_dict[0x25d1f8e46df6eb03]
g2relatives_dict[0xbadd3b54ef02de11]

g2type_dict[0x63af3813d01a9ce3]
g2type_dict[0x528ce72f7e66031b]

g2type_dict[0x5eb06d0946d55b5c]
g2type_dict[0xa7d44ccef97cb142]
g2type_dict[0x46515c67f33148de]

plot_graphene(atom_xy, patoms_collection, indexed_atoms_collection)


g, atom_xy, patoms_collection, indexed_atoms_collection = make_graphene(data_stack[125])
scatter!([g2x_dict[0x25d1f8e46df6eb03],g2x_dict[0xbadd3b54ef02de11]], [g2y_dict[0x25d1f8e46df6eb03], g2y_dict[0xbadd3b54ef02de11]], w=3, xlims=(0,256), ylims=(0,256), aspect_ratio=:equal, leg=false)
plot_atoms(atom_xy)
    scatter(atom_xy[1,:], atom_xy[2,:], w=3, xlims=(0,256), ylims=(0,256), aspect_ratio=:equal, leg=false)



g2p_dict = merge(a2p_dict, b2p_dict, p2p_dict)
g2signature_dict = Dict{UInt, Dict{Int,Int}}()
for (k, v) in g2p_dict
    g2signature_dict[k] = countmap(values(map(x -> g2noa_dict[x], v)))
end

a2p_dict = Dict{UInt, Vector{UInt}}()
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
b2p_dict = Dict{UInt, Vector{UInt}}()
for (k, v) in polygonid_dict
    for vi in last(v)
        if haskey(b2p_dict, bondid_dict_reversed[bond_dict[vi]])
            push!(b2p_dict[bondid_dict_reversed[bond_dict[vi]]],k)
        else
            b2p_dict[bondid_dict_reversed[bond_dict[vi]]] = [k]
        end
    end
end
p2p_dict = Dict{UInt, Vector{UInt}}()
for (k, v) in p2ab_dict
    p2p_dict[k] = unique!(filter!(x -> x != 0 && x != k,collect(Iterators.flatten(map(x -> get(b2p_dict, x, 0), v)))))
end
