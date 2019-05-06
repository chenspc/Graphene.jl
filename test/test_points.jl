using Graphene
using Images
using FileIO
using Plots
using HDF5
using BenchmarkTools

test_file_path = "/Users/chen/Dropbox/_julia/Graphene.jl/test/graphene_defect_nnOutput.h5"
# test_file_path = "/Users/chen/Downloads/graphene_defect_nnOutput.h5"
test_data = h5read(test_file_path, "Dataset1")
# test_image = test_data[:,:,1]
# test_centroids = make_centroids(test_image; threshold=0.5)

# ~300s for 10000 images
@time test_centroids_collection = map(x->make_centroids(test_data[:,:,x]), 1:size(test_data,3))
# ParallelMap
# @time test_centroids_collection = pmap(x->make_centroids(test_data[:,:,x]), 1:size(test_data,3))

centroids2dataframe(test_centroids_collection[1])

plot_centroids(test_centroids_collection[1])
plot_centroids(test_centroids_collection[2])

test_image_path = "/Users/chen/Dropbox/_julia/Graphene.jl/test/test_nnimage.jpeg"
test_gray_im = read_nnimage(test_image_path)
test_centroids = make_centroids(test_gray_im; threshold=0.5)
plot_centroids(test_centroids)

# @benchmark test_centroids = make_centroids(test_gray_im; threshold=0.5)
# @benchmark plot_centroids(test_centroids)






# # Bug? 0 is also considered in component_centroids()
# gray_im = Gray.(load(test_image_path))
# bw_im = gray_im .> 0.5
# markers1 = label_components(erode(bw_im))
# markers2 = label_components(erode(bw_im[100:180,100:180]))
# markers3 = label_components(erode(bw_im[115:165,115:165]))
# markers4 = label_components(erode(bw_im[1:165,1:165]))
# centroids1 = component_centroids(markers1)
# centroids2 = component_centroids(markers2)
# centroids3 = component_centroids(markers3)
# centroids4 = component_centroids(markers4)
# scatter(last.(centroids1[2:end]), first.(centroids1[2:end]), aspect_ratio=:equal, yflip=true, leg=false)
# scatter(last.(centroids2), first.(centroids2), aspect_ratio=:equal, yflip=true, leg=false)
# scatter(last.(centroids3), first.(centroids3), aspect_ratio=:equal, yflip=true, leg=false)
# scatter(last.(centroids4), first.(centroids4), aspect_ratio=:equal, yflip=true, leg=false)
