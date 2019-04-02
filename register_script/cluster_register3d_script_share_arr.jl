using Distributed
if "nodes_name" in keys(ENV)
    nodes_name_str = ENV["nodes_name"]
    node_list = unique(split(nodes_name_str, "\n"))
    machine_spec = [(name, 8) for name in node_list]
    addprocs(machine_spec)
elseif "proc_n" in keys(ENV)
    addprocs(parse(Int,ENV["proc_n"]))
else
    addprocs(8)
end
@everywhere using Printf
print("Process #:", nprocs())
# img_path = "/scratch/binxu.wang/OCPI_data/single_plane_100hz_2min_fish2.imagine"
# save_dir = "/scratch/binxu.wang/OCPI_data/single_plane_100hz_2min_fish2/"
img_path = "/mnt/binxu001/stitched.nhdr"
save_dir = "/mnt/binxu001/stitched/"
if "img_path" in keys(ENV)
    img_path = ENV["img_path"]
end
if "save_dir" in keys(ENV)
    save_dir = ENV["save_dir"]
end
if !isdir(save_dir)
    mkpath(save_dir)
end
@everywhere using LinearAlgebra
@everywhere using Printf, Images
@everywhere mxshift = (10,10,5) #make sure this isn't too small
# mxrot = (0.5,)
# minwidth_rot = fill(0.002, 3)
# SD = Matrix{Float64}(LinearAlgebra.I, 2, 2)
@everywhere begin
    using RegisterQD
    using CoordinateTransformations
    img_path = @fetchfrom(1,img_path)
    img = load(img_path)
    ref_img = img.data[:,:,:,1]
end
@everywhere begin
    function worker(frame_n)
        cur_img = img.data[:,:,:,frame_n]
        @printf("===%d proc processing %d frame=====\n",  myid(), frame_n)
        @time tform, mm = qd_translate(ref_img, cur_img, 
            mxshift; maxevals=500, rtol=0, fvalue=0.0002)
        # qd_rigid(centered(ref_img), centered(cur_img), 
            # mxshift, mxrot, minwidth_rot, SD; maxevals=1000, rtol=0, fvalue=0.0002)
        return tform, mm
    end
end

time_step = size(img,4)
@time result = pmap(worker, 2:10; retry_delays = zeros(5), on_error=ex->(Translation(fill(NaN,3)),NaN)) #time_step
# time_step
# mxrot = (0.01,0.01,0.01)
# minwidth_rot = fill(0.002, 3)
using JLD2, FileIO
using Serialization
array_result = Any[]
for (i, item) in enumerate(result)
    try
        push!(array_result, (item[1].linear, item[1].translation, item[2]))
        if isnan(item[2])
            print("Found NaN result in ", i, "\n")
        end
    catch
        print("Found NaN result in ", i, "\n")
        push!(array_result, NaN)
    end
end
@save joinpath(save_dir, "parallel_reg_result_array.jld2") array_result
serialize(open(joinpath(save_dir, "parallel_reg_result.jls"),"w"),result) # Use IO stream instead of others
# convert the dataformat to more savable one!
