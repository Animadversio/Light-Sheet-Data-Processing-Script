"""
To apply the transformations got from fitting algorithm to really transform the image
(using memory map on disk)
And then do segmentation for the image on disk!
"""
# test_img = testimage("cameraman");
# tfm_arr = deserialize("/Users/binxu/Holy_Optical_Imaging/register_script/parallel_reg_result.jls")
# i = 1
# tfm_img = warp(img,tfm_arr[i][1])
# inds = intersect.(axes(img), axes(tfm_img))
# img[inds...], tfm_img[inds...]
# imshow(img)
orig_img_path = "/mnt/binxu001/03_11_2019/single_plane_100hz_2min_fish2.imagine"
reg_result_path = "/home/binxu/parallel_reg_result.jls"
output_dir = "/mnt/binxu001/03_27_2019/dpf5_fish2_gfp_2d_dynamic/"
if "img_path" in keys(ENV)
    orig_img_path = ENV["img_path"]
end
if "save_dir" in keys(ENV)
    output_dir = ENV["save_dir"]
end
if "reg_result_path" in keys(ENV)
    reg_result_path = ENV["reg_result_path"]
else
    reg_result_path = output_dir * "parallel_reg_result.jls"
end
using Distributed
# if "nodes_name" in keys(ENV)
#     nodes_name_str = ENV["nodes_name"]
#     node_list = unique(split(nodes_name_str, "\n"))
#     machine_spec = [(name, 8) for name in node_list]
#     addprocs(machine_spec)
# else
#     addprocs(6)
# end
addprocs(20)
reg_result_path = "/mnt/binxu001/03_27_2019/dpf5_fish2_gfp_2d_dynamic/parallel_reg_result.jls"
# reg_result_path = "/mnt/binxu001/03_26_2019/40x_4dpf_fish1_gcamp6f/parallel_reg_result.jls"
# orig_img_path = "/mnt/binxu001/03_26_2019/40x_4dpf_fish1_gcamp6f.imagine"
# output_dir = "/mnt/binxu001/03_26_2019/40x_4dpf_fish1_gcamp6f/"
reg_result_path = "/mnt/binxu001/03_26_2019/40x_4dpf_fish1_gcamp6f/parallel_reg_result.jls"
orig_img_path = "/mnt/binxu001/03_26_2019/40x_4dpf_fish1_gcamp6f.imagine"
output_dir = "/mnt/binxu001/03_26_2019/40x_4dpf_fish1_gcamp6f/"
mmap_tmp_file_path = "/mnt/binxu001/tmp_tfm_img2.bin"
orig_img_path = "/mnt/binxu001/stitched.nhdr"
output_dir = "/mnt/binxu001/stitched/"
mmap_tmp_file_path = "/mnt/binxu001/tmp_tfm_img3.bin"
print("Process #:", nprocs())
@everywhere begin
    fit_reg_jls = @fetchfrom 1 reg_result_path
    img_path = @fetchfrom 1 orig_img_path
    tmp_file_path = @fetchfrom 1 mmap_tmp_file_path
    using Images, Printf
    using TestImages, ImageTransformations, CoordinateTransformations, Rotations
    using StaticArrays, Interpolations, LinearAlgebra
    using Serialization #,ImageView
    using SharedArrays
    #using Rebugger
    using FileIO, JLD2
    img = load(img_path)
    tfmed_img_arr = SharedArray{eltype(img),ndims(img)}(tmp_file_path, size(img)
                    ;mode="w+", init=false);
    tfm_arr = deserialize(open(fit_reg_jls, "r"))
    frame_n = size(total_tfm_arr,1)
    if ndims(img) == 3
        total_tfm_arr = pushfirst!([item[1] for item in tfm_arr], AffineMap([1.0 0.0; 0.0 1.0], [0.0, 0.0]))
        inds = axes(centered(img[:,:,1]))
        UL = SVector(inds[1].indices.start, inds[2].indices.start)
        UR = SVector(inds[1].indices.stop, inds[2].indices.start)
        DL = SVector(inds[1].indices.start, inds[2].indices.stop)
        DR = SVector(inds[1].indices.stop, inds[2].indices.stop)
    elseif ndims(img) == 4
        total_tfm_arr = pushfirst!([item[1] for item in tfm_arr], Translation((0.0, 0.0, 0.0)))#IdentityTransformation())
        inds = axes(centered(img[:,:,:,1]))
        UL = SVector(inds[1].indices.start, inds[2].indices.start, inds[3].indices.start)
        UR = SVector(inds[1].indices.stop, inds[2].indices.start, inds[3].indices.start)
        DL = SVector(inds[1].indices.start, inds[2].indices.stop, inds[3].indices.start)
        DR = SVector(inds[1].indices.stop, inds[2].indices.stop, inds[3].indices.start)
        UL_B = SVector(inds[1].indices.start, inds[2].indices.start, inds[3].indices.stop)
        UR_B = SVector(inds[1].indices.stop, inds[2].indices.start, inds[3].indices.stop)
        DL_B = SVector(inds[1].indices.start, inds[2].indices.stop, inds[3].indices.stop)
        DR_B = SVector(inds[1].indices.stop, inds[2].indices.stop, inds[3].indices.stop)
    else
        throw()
    end
    function transform_img(frame_i)
        tfm_img = nothing
        if ndims(img) == 3
            cur_img = img[:,:,frame_i]
            try
                tfm_img =  warp(centered(cur_img), total_tfm_arr[frame_i], inds) # Note we have to transform the image from the center
                print("Transformed frame ", frame_i)
            catch e
                print("Error in registration for frame ",frame_i," ",e)
                tfm_img = cur_img
            end
            tfmed_img_arr[:,:,frame_i] = tfm_img
            # get the 4 corners of the image
            return [total_tfm_arr[frame_i](corner) for corner in [UL, UR, DL, DR]]
        elseif ndims(img) == 4
            cur_img = img[:,:,:,frame_i]
            try
                tfm_img =  warp(centered(cur_img), total_tfm_arr[frame_i], inds) # Note we have to transform the image from the center
                print("Transformed frame ", frame_i)
            catch e
                print("Error in registration for frame ",frame_i," ",e)
                tfm_img = cur_img
            end
            tfmed_img_arr[:,:,:,frame_i] = tfm_img
            # get the 4 corners of the image
            return [total_tfm_arr[frame_i](corner) for corner in [UL, UR, DL, DR, UL_B, UR_B, DL_B, DR_B]]
        end
    end
end
print("Finished defining warp")
@time corner_coord_arr = pmap(transform_img, 1:frame_n) # distributed compute the frame_n
print("Finished transforming ")
serialize(open(joinpath(output_dir, "reg_corner_coord.jls"),"w"),corner_coord_arr)

# all_axes = map(axes,tfmed_img_arr)
# # crop_box = intersect.(all_axes...) # lead to too long processing
# # crop_box = intersect.(map(axes,tfmed_img_arr)...) # lead to stackoverflow!
# crop_box = all_axes[1] # how to find the smallest fitting bbox
# for i = 2:frame_n
#     global crop_box
#     crop_box = intersect.(crop_box, all_axes[i])
# end
# print("Finished getting axes range.")
# @everywhere crop_box_=@fetchfrom 1 crop_box # broadcast the box to every processor
# crop_tfm_img_arr = pmap(frame->frame[crop_box_...], tfmed_img_arr) # distributed compute the crop
# print("Finished cropping image. ")
# img = nothing
# crop_tfm_3d_img = Array{Normed{UInt16,16}}(undef,size(crop_tfm_img_arr[1])...,frame_n)
# for i = 1:frame_n
#     crop_tfm_3d_img[:,:,i] = crop_tfm_img_arr[i]
# end
# print("Finished arranging the cropped image.")
#reshape(collect(crop_tfm_img_arr), (size(crop_tfm_img_arr[1])...,frame_n))
#save("registered_img_arr.jld2","reg_crop_img",crop_tfm_img_arr,"crop_box",crop_box)

# save(joinpath(output_dir, "registered_3dimg.jld2"),"reg_crop_3dimg",tfmed_img_arr,"crop_box",corner_coord_arr)

using CoordinateTransformations
using CellSegmentation, Images, FakeCells, Test, FileIO, JLD2, Statistics
using TileTrees, BoxTrees, IntervalSets
using LinearAlgebra
if ndims(img) == 3
    axis_img = AxisArray(tfmed_img_arr,:y,:x,:time) # assert it's a 2d image
    maxsize = (35,35)
elseif ndims(img) == 4
    axis_img = AxisArray(tfmed_img_arr,:y,:x,:z,:time) # assert it's a 2d image
    maxsize = (35,35,5)
end
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
