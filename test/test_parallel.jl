using Distributed
addprocs(4)
nworkers()
nprocs()

@everywhere using Pkg; Pkg.activate("/Users/chen/.julia/environments/v1.1/")
@everywhere using Graphene, FileIO, HDF5, JLD, JuliaDB, CSV, Test, BenchmarkTools, Profile, Images
@everywhere using SharedArrays
using Test
using BenchmarkTools

test_mat = ones(5,5,5)

test_new_mat = [test_mat[:,:,i] .+ i for i in 1:5]

@benchmark test_new_mat = pmap(x -> x^20, 1:500)
@benchmark test_new_mat = map(x -> x^20000, 1:500)

@time test_new_mat = map(x -> x^20000, 1:500^3)

@distributed for N in 1:20
           println("The N of this iteration in $N")
           @time test_new_mat = map(x -> x^20000, 1:500^3)
       end

test_dict = Dict{Int, Int}()
test_dict_new = Dict{Int, Int}()

for i in 1:10000
    test_dict[i] = 10000 + i
end

test_dict
length(test_dict)

@distributed merge for i in 1:length(test_dict)
    test_dict_new[i] = test_dict[i] + 100000
end

test_dict
test_dict_new

ntails = @distributed (+) for i=1:2000000
   rand(Bool)
end

ntails

pmap(i in 1:length(test_dict)
    test_dict_new[i] = test_dict[i] + 100000
end


M = Matrix{Float64}[rand(1000,1000) for i = 1:10]
pmap(sum, M)

image_series = Array{Float64}[ones(1000,1000) for i = 1:10]
pmap(sum, image_series)
image_series_dict = Dict([i => ones(1000,1000) for i = 1:10])
pmap(sum, image_series_dict)


a = SharedArray{Float64}(10)
# b = SharedArray{Float64}(rand(5,5,10))
b = rand(5,5,10)
# c = Array{Float64}[ones(5,5) for i = 1:10]
c = [ones(6,6)*i for i = 1:10]
@distributed for i = 1:10
    # a[i] = sum(b[:,:,i])
    @time a[i] = sum(c[i])*1.00004523452001^1024541
    println("Result: $(a[i]) The i of this iteration in $i")
end

@show a

Sys.CPU_THREADS





@time data_stack = import_stack("/Users/chen/Dropbox/_julia/Graphene.jl/test/graphene_defect_nnOutput_1.h5")

@everywhere data_stack_small = $data_stack[1:100]
a = pmap(x -> make_graphene(data_stack_small[x]; frame=x), collect(1:length(data_stack_small)))

@time g_stack = map(x -> make_graphene(data_stack[x]; frame=x), collect(1:length(data_stack)))
@time g_stack = pmap(x -> make_graphene(data_stack[x]; frame=x), [1:100])
@time g_stack = map(x -> make_graphene(data_stack[x]; frame=x), collect(1:length(data_stack[1:100])))

data_stack[1:100]
data_stack
