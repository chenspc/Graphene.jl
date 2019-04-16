using DataFrames
using DataFramesMeta

using Images
using ImageMorphology
using ImageSegmentation
using ImageView
using TestImages
# using QuartzImageIO
using ImageMagick
using FileIO
using Plots

export read_nnimage, make_centroids, plot_centroids

function read_nnimage(image_path)
    gray_im = Gray.(load(image_path))
end

function make_centroids(gray_im; threshold=0.5)
    bw_im = gray_im .> threshold
    markers = label_components(erode(bw_im))
    centroids = component_centroids(markers)
end

function plot_centroids(centroids)
    scatter(last.(centroids), first.(centroids), aspect_ratio=:equal, yflip=true, leg=false)
end
