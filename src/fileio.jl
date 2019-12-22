export import_csv
export read_nnimage
export make_centroids
export centroid2xy
export import_stack
export stack2xy

function import_csv(centroid_data_file::String)
    nnResult = read(centroid_data_file; header=false)
end

function read_nnimage(image_path::String)
    gray_im = Gray.(load(image_path))
end

function make_centroids(gray_im; threshold=0.5)
    bw_im = gray_im .> threshold
    markers = label_components(erode(bw_im))
    centroids = component_centroids(markers)
    popfirst!(centroids)
    return centroids
end

function centroid2xy(centroids)
    xy = [first.(centroids) last.(centroids)]
    atom_xy = permutedims(xy, (2,1))
end

function import_stack(nnResult_stack_path::String; range=[])
    if range == []
        nnResult_stack = h5read(nnResult_stack_path, "Dataset1")
    else
        nnResult_stack = h5read(nnResult_stack_path, "Dataset1", (:,:,range))
    end
    return nnResult_stack
end

function stack2xy(nnResult_stack)
    centroids_collection = map(x->make_centroids(nnResult_stack[:,:,x]), 1:size(nnResult_stack,3))
    xy_stack = map(centroid2xy, centroids_collection)
    return xy_stack
end
