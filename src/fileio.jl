export import_csv
export data_reshape
export import_stack
export centroids2xy

function import_csv(centroid_data_file::String)
    nnResult = CSV.read(centroid_data_file; header=false)
end

function data_reshape(atom_data; image_sampling=1)
    atom_xy = convert(Array{Float64}, Matrix(atom_data)') .* image_sampling
end

function import_stack(nnResult_stack_path::String; range=[])
    if range == []
        nnResult_stack = h5read(nnResult_stack_path, "Dataset1")
    else
        nnResult_stack = h5read(nnResult_stack_path, "Dataset1", (:,:,range))
    end
    centroids_collection = map(x->make_centroids(nnResult_stack[:,:,x]), 1:size(nnResult_stack,3))
    data_stack = map(centroids2xy, centroids_collection)
    return data_stack
end

function centroids2xy(centroids)
    xy = [first.(centroids) last.(centroids)]
    atom_xy = permutedims(xy, (2,1))
end
