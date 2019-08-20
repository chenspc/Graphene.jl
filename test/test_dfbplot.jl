using CSV
using Images
using Gadfly

dfb_data = CSV.read("/Users/chen/Dropbox/_julia/Graphene.jl/test/processing/20180517 121632/graphene_defect_nnOutput0003_dfb_analysis.csv")

dfb_data
d_data = dfb_data[dfb_data[:, :noa] .== 14, :]
f_data = dfb_data[dfb_data[:, :noa] .== 1, :]
b_data = dfb_data[dfb_data[:, :noa] .== 6, :]
display_dfb(d_data, f_data, b_data; shape="line")
display_dfb(d_data, f_data, b_data)
