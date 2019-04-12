export import_csv
export data_reshape

function import_csv(centroid_data_file::String)
    raw_data = CSV.read(centroid_data_file; header=false)
end

function data_reshape(atom_data; image_sampling=1)
    atom_xy = convert(Array{Float64}, Matrix(atom_data)') .* image_sampling
end
