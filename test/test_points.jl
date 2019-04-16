using Graphene
using Images
using ImageMorphology
using ImageSegmentation
using ImageView
using TestImages
# using QuartzImageIO
using ImageMagick
using FileIO
using Plots

test_image_path = "/Users/chen/Dropbox/_julia/Graphene.jl/test/test_nnimage.jpeg"

markers = label_components(ImageMorphology.erode(test_bw_image))

segments = watershed(Gray.(test_image), markers)
segments = watershed(test_bw_image, markers)
atom_positions = component_centroids(markers)

ImageMorphology.opening(Gray.(test_bw_image))
ImageMorphology.erode(Gray.(test_bw_image))
ImageSegmentation.watershed(test_bw_image)

ImageView.closeall()


test_image_path = "/Users/chen/Dropbox/_julia/Graphene.jl/test/test_nnimage.jpeg"
test_gray_im = read_nnimage(test_image_path)
test_gray_im = test_gray_im[100:200,100:200]
test_centroids = make_centroids(test_gray_im; threshold=0.5)
plot_centroids(test_centroids)

test_bw_im = test_gray_im .> 0.5

Gray.(test_bw_im[115:165, 115:165])

Gray.(test_gray_im[115:165, 115:165])

Gray.(erode(test_bw_im)[115:165, 115:165])

Gray.(erode(test_bw_im))

test_markers = label_components(Gray.(erode(test_bw_im)))
test_markers[115:165, 115:165]
test_markers[115:165, 115:165]
Gray.(label_components(Gray.(erode(test_bw_im[115:165, 115:165]))) .> 0)

component_centroids(test_markers)
plot_centroids(component_centroids(test_markers))

Gray.(watershed(opening(test_bw_im)))
