using CoordinateTransformations
using CellSegmentation, Images, FakeCells, Test, FileIO, JLD2, Statistics
using TileTrees, BoxTrees, IntervalSets
using LinearAlgebra
img_path = "/scratch/binxu.wang/OCPI_data/single_plane_100hz_2min_fish2.imagine"
#"/mnt/binxu001/registered_3dimg.jld2"
output_dir = "/scratch/binxu.wang/OCPI_data/single_plane_100hz_2min_fish2"
#"/mnt/binxu001/03_11_2019/proc_100hz_2min_fish2"
if "save_dir" in keys(ENV)
    output_dir = ENV["save_dir"]
end
if "reg_img_path" in keys(ENV)
    img_path = ENV["reg_img_path"]
else
    img_path = output_dir * "registered_3dimg.jld2"
end
if occursin(".jld2", img_path)
    f_dict = load(img_path)
    img = f_dict["reg_crop_3dimg"]
elseif occursin(".imagine", img_path)
    img = load(img_path)
else
    throw(ErrorException)
end
if !isdir(output_dir)
    mkpath(output_dir)
end

axis_img = AxisArray(img,:y,:x,:time)
maxsize = (35,35)
tempfile = tempname() * ".jld2" # "/scratch/binxu.wang/temp_tree_data1.jld2"#
@time segment_prepare(tempfile, axis_img, maxsize, 0) # ; ncomp=16
# tempfile = "/scratch/binxu.wang/OCPI_data/temp_tree_data.jld2"
@time ttree, reftree = segment_merge(tempfile)
rm(tempfile)
using Serialization
serialize(open(joinpath(output_dir,"segment_ttree.jls"),"w"),ttree)
serialize(open(joinpath(output_dir,"segment_reftree.jls"),"w"),reftree)
save(joinpath(output_dir,"segment_tree.jld2"),"ttree",ttree,"reftree",reftree)
#%%
# load("/scratch/binxu.wang/segment_tree.jld2","ttree")
# load("/scratch/binxu.wang/temp_tree_data.jld2")
