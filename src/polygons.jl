export make_polygons
export path2turn
export link_turns!

function make_polygons(paths)
    turn_dict = Dict(map(path2turn, collect(paths)))
    collector = []
    while !isempty(turn_dict)
        push!(collector, link_turns!(turn_dict))
    end
    return collector
end

function path2turn(path)
    (path[1], path[2]) => (path[2], path[3])
end

function link_turns!(turn_dict)
    carrier = pop!(turn_dict, first(first(turn_dict)), "empty")
    polygon_atoms = [first(carrier)]
    polygon_bonds = [carrier]
    while (carrier = pop!(turn_dict, carrier, "empty")) != "empty"
        push!(polygon_atoms, first(carrier))
        push!(polygon_bonds, carrier)
    end
    return polygon_atoms, polygon_bonds
end
