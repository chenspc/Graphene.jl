
export id_generator

function id_generator(c::Channel)
# All-in-one ID generator for AtomID, BondID, and PolygonID
    alphabet = collect('a':'z')
    for n1 in alphabet
        for n2 in alphabet
            for n3 in alphabet
                for n4 in alphabet
                    put!(c, string(n1,n2,n3,n4))
                end
            end
        end
    end
end
