using Graphene
using FileIO
using HDF5
using JLD
using CSV
using Test
using BenchmarkTools
using Profile
# using ProfileView

data = import_csv("/Users/chen/Dropbox/_julia/Graphene.jl/test/test_data_flower.csv")
@benchmark g = make_graphene(data)
g = make_graphene(data)

@time data_stack = import_stack("/Users/chen/Dropbox/_julia/Graphene.jl/test/graphene_defect_nnOutput.h5")
@time g_stack = map((i, x) -> make_graphene(x;frame=i, dataset="test_dataset"), collect(1:300),data_stack[1:300])
g_merge = g_stack[1]
map(x -> merge!(x, g_merge; pkey=:frame), g_stack[2:300])
g_merge
merge(g_stack[3], [])

length(g_stack)
function merge_tables(t_stack)
    t_merge = t_stack[1]
    for i = 2:length(t_stack)
        t_merge = merge(t_stack[i], t_merge)    
    end
end
@time t_merge = merge_tables(g_stack)

# 3858.473570 seconds (2.93 G allocations: 1.290 TiB, 72.25% gc time)
make_graphene(data_stack[125])
make_graphene(data_stack[80])
make_graphene(data_stack[81])
make_graphene(data_stack[82])

save("/Users/chen/Downloads/test_gaResult.jld", "g_stack", g_stack)
save("/Users/chen/Downloads/test_gaResult.jld", "t", t, "arr", z)

save(File(format"JLD","/Users/chen/Downloads/test_gaResult.jld"), "g_stack", g_stack)
save(File(format"JLD","/Users/chen/Downloads/test_gaResult_Tuple.jld"), "g_stack", g_stack)

t = 15
z = [1,3]
save("/tmp/myfile2.jld", "t", t, "arr", z)

@profile make_graphene(data_stack[1])
Profile.print()
ProfileView.view()

save("/Users/chen/Downloads/test_gaResult_new.csv", g)
save("/Users/chen/Downloads/test_gaResult_new_stack.csv", merge(g_stack[1],g_stack[2]))
loadtable("/Users/chen/Downloads/test_gaResult.csv")
a = CSV.read("/Users/chen/Downloads/test_gaResult_new.csv")
