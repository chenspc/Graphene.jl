using Graphene
using Glob
using CSV
using Plots

a_list = glob("*/*.csv", "/Users/chen/Dropbox/_julia/Graphene.jl/test/processing/")

for i in a_list
    a = CSV.read(i)
    d = a[findall(in([14, 8]), a.noa), :]
    f = a[findall(in([1]), a.noa), :]
    b = a[findall(in([6]), a.noa), :]
   display_dfb(d, f, b)
end

for i in a_list
    i
end
display_dfb(d, f, b)

plotkkkkk
