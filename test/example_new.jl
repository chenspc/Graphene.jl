using Graphene
using FileIO
using HDF5
using JLD
using JuliaDB
using CSV
using Test
using BenchmarkTools
using Profile
using Images
# using ProfileView
using MAT

data = import_csv("/Users/chen/Dropbox/_julia/Graphene.jl/test/test_data_flower.csv")
@benchmark g = make_graphene(data)
g = make_graphene(data)
# FileIO.save("/Users/chen/Downloads/test_gaResult.csv", g)

@time data_stack = import_stack("/Users/chen/Dropbox/_julia/Graphene.jl/test/graphene_defect_nnOutput.h5")
# @time g_stack = map(x -> make_graphene(x), data_stack[1:300])
# @time g_stack = map(x -> make_graphene(data_stack[x]; frame=x), collect(4801:4850))
# @time g_stack_bug = map(x -> make_graphene(data_stack[x]; frame=x), collect(29:30))

@time g_stack = map(x -> make_graphene(data_stack[x]; frame=x), collect(1:length(data_stack)))
# 2280.918385 seconds (3.05 G allocations: 1.326 TiB, 48.70% gc time) 20190514
@time d_stack = map(x -> find_defect(x, "all"), g_stack)
# 9.768841 seconds (85.64 M allocations: 3.683 GiB, 43.98% gc time)
@time d_flower_stack = map(x -> find_defect(x, "flower"), g_stack)
# 6.680508 seconds (84.68 M allocations: 3.618 GiB, 33.14% gc time)
@time d_butterfly_stack = map(x -> find_defect(x, "butterfly"), g_stack)
# 6.721757 seconds (84.69 M allocations: 3.619 GiB, 33.37% gc time)
@time d_divacancy_stack = map(x -> find_defect(x, "divacancy"), g_stack)
# 7.000525 seconds (84.68 M allocations: 3.618 GiB, 34.73% gc time)

@time d_merge = merge_stack(d_stack)
# 12.554718 seconds (93.69 M allocations: 4.699 GiB, 21.72% gc time)
@time d_flower_merge = merge_stack(d_flower_stack)
# 0.905638 seconds (3.39 M allocations: 443.802 MiB, 24.71% gc time)
# "AAAAAAAAAwAAAAAAAAAAAAAAAAAA"
@time d_butterfly_merge = merge_stack(d_butterfly_stack)
# 7.409825 seconds (58.27 M allocations: 2.970 GiB, 22.43% gc time)
# "AAAAAAIABAAAAAAAAAAAAAAAAAAA"
@time d_divacancy_merge = merge_stack(d_divacancy_stack)
# 0.613168 seconds (1.99 M allocations: 274.278 MiB, 21.33% gc time)
# "AAAAAAAKAAAAAAAAAAAAAAAAAAAA"

d_flower_stack_isolated = map(x -> find_isolated_defect(x, "flower"), g_stack)
d_butterfly_stack_isolated = map(x -> find_isolated_defect(x, "butterfly"), g_stack)
d_divacancy_stack_isolated = map(x -> find_isolated_defect(x, "divacancy"), g_stack)

d_flower_merge_isolated = merge_stack(d_flower_stack_isolated)
d_butterfly_merge_isolated = merge_stack(d_butterfly_stack_isolated)
d_divacancy_merge_isolated = merge_stack(d_divacancy_stack_isolated)

display_dfb(d_divacancy_merge, d_flower_merge, d_butterfly_merge)
display_dfb(d_divacancy_merge, d_flower_merge, d_butterfly_merge; shape="line")
display_dfb(d_divacancy_merge_isolated, d_flower_merge_isolated, d_butterfly_merge_isolated)
display_dfb(d_divacancy_merge_isolated, d_flower_merge_isolated, d_butterfly_merge_isolated; shape="line")

a = display_dfb(d_divacancy_merge_isolated, d_flower_merge_isolated, d_butterfly_merge_isolated)

FileIO.save("~/Downloads/test.png", display_dfb(d_divacancy_merge_isolated, d_flower_merge_isolated, d_butterfly_merge_isolated))

FileIO.save("/Users/chen/Dropbox/_julia/Graphene.jl/examples/d_merge.csv", d_merge)
FileIO.save("/Users/chen/Dropbox/_julia/Graphene.jl/examples/d_flower_merge.csv", d_flower_merge)
FileIO.save("/Users/chen/Dropbox/_julia/Graphene.jl/examples/d_butterfly_merge.csv", d_butterfly_merge)
FileIO.save("/Users/chen/Dropbox/_julia/Graphene.jl/examples/d_divacancy_merge.csv", d_divacancy_merge)

# FileIO.save("/Users/chen/Dropbox/_julia/Graphene.jl/examples/d_merge_isolated.csv", d_merge_isolated)
FileIO.save("/Users/chen/Dropbox/_julia/Graphene.jl/examples/d_flower_merge_isolated.csv", d_flower_merge_isolated)
FileIO.save("/Users/chen/Dropbox/_julia/Graphene.jl/examples/d_butterfly_merge_isolated.csv", d_butterfly_merge_isolated)

# @time GFrame
# 3858.473570 seconds (2.93 G allocations: 1.290 TiB, 72.25% gc time)
# make_graphene(data_stack[125])
# make_graphene(data_stack[80])
# make_graphene(data_stack[81])
# make_graphene(data_stack[82])
loadtable("/Users/chen/Downloads/test_gaResult.csv")
a = CSV.read("/Users/chen/Downloads/test_gaResult.csv")

# @profile make_graphene(data_stack[1])
# Profile.print()
# ProfileView.view()
#
@everywhere using Pkg; Pkg.activate(".")
@everywhere using Graphene, FileIO, HDF5, JLD, JuliaDB, CSV, Test, BenchmarkTools, Profile, Images
@everywhere data_stack_small = $data_stack[1:100]
pmap(x -> make_graphene(data_stack_small[x]; frame=x), collect(1:length(data_stack_small)))
