module Graphene


# All-in-one ID generator for AtomID, BondID, and PolygonID
function id_generator(c::Channel)
    for n1 in collect('a':'z')
        for n2 in collect('a':'z')
            for n3 in collect('a':'z')
                for n4 in collect('a':'z')
                    put!(c, string(n1,n2,n3,n4))
                end
            end
        end
    end
end

end # module
