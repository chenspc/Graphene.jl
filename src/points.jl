# using DataFrames
# using DataFramesMeta
# using Images
# using FileIO
# using Plots

export read_nnimage, make_centroids, plot_centroids, centroids2dataframe

function read_nnimage(image_path)
    gray_im = Gray.(load(image_path))
end

function make_centroids(gray_im; threshold=0.5)
    bw_im = gray_im .> threshold
    markers = label_components(erode(bw_im))
# component_centroids() is currently the bottleneck for the whole workflow
    centroids = component_centroids(markers)
    popfirst!(centroids)
    return centroids
end

function plot_centroids(centroids)
    scatter(last.(centroids), first.(centroids), aspect_ratio=:equal, yflip=true, leg=false)
end

function centroids2dataframe(centroids)
    df = DataFrame(Column1 = first.(centroids), Column2 = last.(centroids))
    return df
end
