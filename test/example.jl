using Graphene
using FileIO
using HDF5
using JLD
using Test
using BenchmarkTools
using Profile
using ProfileView

data = import_csv("/Users/chen/Dropbox/_julia/Graphene.jl/test/test_data_flower.csv")
@benchmark g = make_graphene(data)

@time data_stack = import_stack("/Users/chen/Dropbox/_julia/Graphene.jl/test/graphene_defect_nnOutput.h5")
@time g_stack = map(x -> make_graphene(x), data_stack[1:300])
# 3858.473570 seconds (2.93 G allocations: 1.290 TiB, 72.25% gc time)
make_graphene(data_stack[125])
make_graphene(data_stack[80])
make_graphene(data_stack[81])
make_graphene(data_stack[82])

save("/Users/chen/Downloads/test_gaResult.jld", "g_stack", g_stack)
save("/Users/chen/Downloads/test_gaResult.jld", "t", t, "arr", z)

save(File(format"JLD","/Users/chen/Downloads/test_gaResult.jld"), "g_stack", g_stack)

t = 15
z = [1,3]
save("/tmp/myfile2.jld", "t", t, "arr", z)

@profile make_graphene(data_stack[1])
Profile.print()
ProfileView.view()
