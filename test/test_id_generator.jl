using Test
using Graphene

@test 1 == 1
@test_broken 1 == 0

chnl_test = Channel(id_generator)
id_test = []
for i in 1:30
        push!(id_test, take!(chnl_test))
end
@test id_test == ["aaaa", "aaab", "aaac", "aaad", "aaae", "aaaf", "aaag", "aaah", "aaai", "aaaj",
                  "aaak", "aaal", "aaam", "aaan", "aaao", "aaap", "aaaq", "aaar", "aaas", "aaat",
                  "aaau", "aaav", "aaaw", "aaax", "aaay", "aaaz", "aaba", "aabb", "aabc", "aabd"]
