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
img_path = "/scratch/binxu.wang/OCPI_data/single_plane_100hz_2min_fish2.imagine"
save_dir = "/scratch/binxu.wang/OCPI_data/single_plane_100hz_2min_fish2/"
if "img_path" in keys(ENV)
    img_path = ENV["img_path"]
end
if "save_dir" in keys(ENV)
    save_dir = ENV["save_dir"]
end
if !isdir(save_dir)
    mkpath(save_dir)
end

using Images
img = load(img_path)
@everywhere begin
    using RegisterQD
    using Images
    using CoordinateTransformations
    using LinearAlgebra
    # using SharedArrays
    using Printf
    SD = Matrix{Float64}(LinearAlgebra.I, 2, 2)
    mxshift = (100,100) #make sure this isn't too small
    mxrot = (0.5,)
    minwidth_rot = fill(0.002, 3)
    ref_img = @fetchfrom 1 img[:,:,1]
    function worker(frame_n)
        cur_img = @fetchfrom 1 img[:,:,frame_n]
        @printf("===%d proc processing %d frame=====",  myid(), frame_n)
        @time tform, mm =qd_rigid(centered(ref_img), centered(cur_img), 
            mxshift, mxrot, minwidth_rot, SD; maxevals=1000, rtol=0, fvalue=0.0002)
        return tform, mm
    end
end

time_step = size(img,3)
@time result = pmap(worker, 2:time_step; retry_delays = zeros(5), on_error=ex->(AffineMap(fill(NaN, (2,2)), fill(NaN, 2)),NaN)) #time_step

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
