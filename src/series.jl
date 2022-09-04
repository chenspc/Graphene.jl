export make_gsa
export link_past_future!
export chain_time_series!
export plot_graphene_gsa
export save2gif

"""
    make_gsa(graphene::AbstractVector{<:AbstractGEntry})

Make a StructArray from a Graphene object, `graphene`, for efficient computation.
"""
function make_gsa(graphene::AbstractVector{<:AbstractGEntry})
    sa_id        = get_id.(graphene)
    sa_x         = get_x.(graphene)
    sa_y         = get_y.(graphene)
    sa_relatives = get_relatives.(graphene)
    sa_signature = get_signature.(graphene)
    sa_frame     = get_frame.(graphene)
    sa_dataset   = get_dataset.(graphene)
    sa_noa       = get_noa.(graphene)
    sa_members   = get_members.(graphene)
    sa_type      = get_type.(graphene)
    sa_past      = get_past.(graphene)
    sa_future    = get_future.(graphene)

    graphene_sa = StructArray(id        = sa_id,
                              x         = sa_x,
                              y         = sa_y,
                              relatives = sa_relatives,
                              signature = sa_signature,
                              frame     = sa_frame,
                              dataset   = sa_dataset,
                              noa       = sa_noa,
                              members   = sa_members,
                              type      = sa_type,
                              past      = sa_past,
                              future    = sa_future)
    return graphene_sa
end

"""
    link_past_future!(gpast::AbstractVector, gfuture::AbstractVector; Δxy=5, atoms_only=false)

Link the atoms, bonds, polygons and defects in two Graphene objects, `gpast` and `gfuture`. `Δxy` is the maximum distance for two items to be considered related. 

If `atoms_only` is set to `true`, the code will only link the atoms but not for all other types of items. 

If one atom in `gpast` has two or more possible corresponding atoms in `gfuture`, it will link the closest pair.

The pairing process is recursive. 

Internally, the future or past of an item is represented as a complex number, where the real part is the ID and the imaginary part is the frame number.  
"""
function link_past_future!(gpast::AbstractVector, gfuture::AbstractVector; Δxy=5, atoms_only=false)
    gpast_atoms = view(gpast, [i for i in 1:length(gpast) if gpast.type[i] == :atom && gpast.future[i] == 0im])
    gfuture_atoms = view(gfuture, [i for i in 1:length(gfuture) if gfuture.type[i] == :atom && gfuture.past[i] == 0im])

    if !isempty(gpast_atoms) && !isempty(gfuture_atoms)
        gpast_xy = hcat(gpast_atoms.x, gpast_atoms.y)'
        gfuture_xy = hcat(gfuture_atoms.x, gfuture_atoms.y)'

        kdtree_gpast = KDTree(gpast_xy, leafsize = 10)
        kdtree_gfuture = KDTree(gfuture_xy , leafsize = 10)

        last_round = true

        possible_future, future_shift = knn(kdtree_gfuture, gpast_xy, 1, true)
        possible_past, past_shift = knn(kdtree_gpast, gfuture_xy, 1, true)

        for i in 1:length(possible_future)
            if  !isempty(possible_future[i]) && first(future_shift[i]) ≤ Δxy && i == first(possible_past[first(possible_future[i])])

                gpast_atoms.future[i] = first(possible_future[i]) |> x -> gfuture_atoms.id[x] + gfuture_atoms.frame[x]*im
                gfuture_atoms.past[first(possible_future[i])] = i |> x -> gpast_atoms.id[x] + gpast_atoms.frame[x]*im

                last_round = false
            end
        end

        last_round || link_past_future!(gpast, gfuture; Δxy=Δxy, atoms_only=true)
    end

    if atoms_only == false
        gpast_others = view(gpast, [i for i in 1:length(gpast) if gpast.type[i] != :atom])
        map(gpast_others) do x
            atom_relatives= view(gpast, [i for i in x.relatives if gpast.type[i] == :atom])
            future_atom_ids = real(atom_relatives.future)
            if !|(iszero.(future_atom_ids)...)
                atom_relatives_future = view(gfuture, future_atom_ids)
                possible_future_other = view(gfuture, [i for i in intersect(atom_relatives_future.relatives...) if gfuture.type[i] == x.type])
                for p in 1:length(possible_future_other)
                    possible_future_other.past[p] == 0im || continue
                    gpast.future[x.id] = possible_future_other.id[p] + possible_future_other.frame[p]*im
                    possible_future_other.past[p] = x.id + x.frame*im
                    break
                end
            end
        end
    end
    return nothing
end

function chain_time_series!(gseries::AbstractVector; Δt_max=1, kwargs...)
    for i in 1:Δt_max
        map((p,f) -> link_past_future!(p, f; kwargs...), view(gseries, 1:length(gseries)-i), view(gseries, 1+i:length(gseries)))
    end
    return nothing
end

function plot_graphene_gsa(gsa)
    gsa_atoms  = view(gsa, [i for i in 1:length(gsa) if gsa.type[i] == :atom])
    links      = view(gsa_atoms, [i for i in 1:length(gsa_atoms) if gsa_atoms.past[i] != 0im && gsa_atoms.future[i] != 0im])
    heads      = view(gsa_atoms, [i for i in 1:length(gsa_atoms) if gsa_atoms.past[i] == 0im && gsa_atoms.future[i] != 0im])
    tails      = view(gsa_atoms, [i for i in 1:length(gsa_atoms) if gsa_atoms.past[i] != 0im && gsa_atoms.future[i] == 0im])
    singletons = view(gsa_atoms, [i for i in 1:length(gsa_atoms) if gsa_atoms.past[i] == 0im && gsa_atoms.future[i] == 0im])
    scatter(links.x, links.y, aspect_ratio=1, label="links")
    scatter!(heads.x, heads.y, aspect_ratio=1, label="heads")
    scatter!(tails.x, tails.y, aspect_ratio=1, label="tails")
    scatter!(singletons.x, singletons.y, aspect_ratio=1, label="singletons")
end

function save2gif(gif_path, gsa; fps=10)
    anim = @animate for i in gsa
        plot_graphene_gsa(i)
    end
    gif(anim, gif_path, fps=fps)
end
