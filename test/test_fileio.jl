using Graphene
using Test
using DataFrames
using Images: Gray
using FixedPointNumbers
using HDF5: h5read

test_csv_data = import_csv("test_files/test_data_flower.csv")
@test isa(test_csv_data, DataFrame)
@test ncol(test_csv_data) == 2

@test isa(read_nnimage("test_files/test_nnimage_flower.jpeg"), Array{Gray{FixedPointNumbers.Normed{UInt8,8}},2})

test_nnResult_stack = import_stack("test_files/graphene_defect_nnOutput.h5", range=1:10)
@test size(test_nnResult_stack) == (256, 256, 10)
test_nnResult_stack = import_stack("test_files/graphene_defect_nnOutput.h5")
@test size(test_nnResult_stack) == (256, 256, 100)

test_centroids = make_centroids(test_nnResult_stack[:,:,1])
@test isa(test_centroids, Array{Tuple{Float64,Float64},1})

@test isa(centroid2xy(test_centroids), Array{Float64,2})
@test size(centroid2xy(test_centroids)) == (2, 837)

@test isa(stack2xy(test_nnResult_stack), Array{Array{Float64,2},1})
