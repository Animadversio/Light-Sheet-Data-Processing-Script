"""
To apply the transformations got from fitting algorithm to really transform the image
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
output_dir = "/mnt/binxu001/03_11_2019/proc_100hz_2min_fish2"
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
addprocs(7)
using Images, Printf
img = load(orig_img_path)
print("Process #:", nprocs())
@everywhere begin
    fit_reg_jls = @fetchfrom 1 reg_result_path
    using TestImages, Images, ImageTransformations, CoordinateTransformations, Rotations
    using StaticArrays, Interpolations, LinearAlgebra
    using Serialization #,ImageView
    #using Rebugger
    using FileIO, JLD2
    tfm_arr = deserialize(open(fit_reg_jls, "r"))
    total_tfm_arr = pushfirst!([item[1] for item in tfm_arr], AffineMap([1.0 0.0; 0.0 1.0], [0.0, 0.0]))#IdentityTransformation())
    frame_n = size(total_tfm_arr,1)
    function transform_img(frame_i)
        cur_img = @fetchfrom 1 img[:,:,frame_i]
        try
            return warp(centered(cur_img), total_tfm_arr[frame_i]) # Note we have to transform the image from the center
        catch e
            print("Error in registration for frame ",frame_i," ",e)
            return cur_img
        end
    end
end
print("Finished defining warp")
tfmed_img_arr = pmap(transform_img, 1:frame_n) # distributed compute the
print("Finished transforming ")
#[ warp(img[:,:,i], total_tfm_arr[i]) for i in 1:frame_n]
all_axes = map(axes,tfmed_img_arr)
# crop_box = intersect.(all_axes...) # lead to too long processing
# crop_box = intersect.(map(axes,tfmed_img_arr)...) # lead to stackoverflow!
crop_box = all_axes[1] # how to find the smallest fitting bbox
for i = 2:frame_n
    global crop_box
    crop_box = intersect.(crop_box, all_axes[i])
end
print("Finished getting axes range. ")
@everywhere crop_box_=@fetchfrom 1 crop_box # broadcast the box to every processor
crop_tfm_img_arr = pmap(frame->frame[crop_box_...], tfmed_img_arr) # distributed compute the crop
print("Finished cropping image. ")
img = nothing
crop_tfm_3d_img = Array{Normed{UInt16,16}}(undef,size(crop_tfm_img_arr[1])...,frame_n)
for i = 1:frame_n
    crop_tfm_3d_img[:,:,i] = crop_tfm_img_arr[i]
end
print("Finished arranging the cropped image.")
#reshape(collect(crop_tfm_img_arr), (size(crop_tfm_img_arr[1])...,frame_n))
#save("registered_img_arr.jld2","reg_crop_img",crop_tfm_img_arr,"crop_box",crop_box)
save(joinpath(output_dir, "registered_3dimg.jld2"),"reg_crop_3dimg",crop_tfm_3d_img,"crop_box",crop_box)
