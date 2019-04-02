using CoordinateTransformations
using CellSegmentation, Images, FakeCells, Test, FileIO, JLD2, Statistics
using TileTrees, BoxTrees, IntervalSets
using LinearAlgebra

# img_path = "/mnt/binxu001/03_11_2019/single_plane_100hz_2min_fish2.imagine"
img_path = "/mnt/binxu001/registered_3dimg.jld2"
output_dir = "/mnt/binxu001/03_11_2019/proc_100hz_2min_fish2"
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
axis_img = AxisArray(img,:y,:x,:time) # assert it's a 2d image
maxsize = (35,35)
tempfile = tempname() * ".jld2" # "temp_tree_data.jld2"
@time segment_prepare(tempfile, axis_img, maxsize, 0) # ; ncomp=16
@time ttree, reftree = segment_merge(tempfile)
rm(tempfile)
using Serialization
serialize(open(joinpath(output_dir,"segment_ttree.jls"),"w"),ttree)
serialize(open(joinpath(output_dir,"segment_reftree.jls"),"w"),reftree)
# serialize(open("/scratch/binxu.wang/OCPI_data/segment_ttree.jls","w"),ttree)
# serialize(open("/scratch/binxu.wang/OCPI_data/segment_reftree.jls","w"),reftree)
save(joinpath(output_dir,"segment_tree.jld2"),"ttree",ttree,"reftree",reftree)
# load("/mnt/binxu001/segment_tree.jld2","ttree","reftree")
# load("temp_tree_data.jld2")

"""
for the (779, 797, 12000) image shape 
Running speed
segment_prepare:
[ Info: Saving to disk, please wait
1122.210274 seconds (317.07 M allocations: 299.542 GiB, 4.36% gc time)
segment_merge:
After centering: 449
251.564436 seconds (79.74 M allocations: 56.772 GiB, 7.06% gc time)
"""
