using Distributed
#addprocs(10)
# addprocs([("node006",8),("node007",8),("node008",8),("node009",8),("node010",8)])
nodes_name_str = ENV["nodes_name"]
node_list = unique(split(nodes_name_str, "\n"))
machine_spec = [(name, 8) for name in node_list]
addprocs(machine_spec)
@everywhere using Printf
print("Process #:", nprocs())
# using Images
# img = load("/scratch/binxu.wang/OCPI_data/single_plane_100hz_2min_fish2.imagine")

@everywhere begin
    using RegisterQD
    using Images
    using CoordinateTransformations
    using LinearAlgebra
    # using SharedArrays
    using Printf
    img = load("/scratch/binxu.wang/OCPI_data/single_plane_100hz_2min_fish2.imagine")
    SD = Matrix{Float64}(LinearAlgebra.I, 2, 2)
    mxshift = (100,100) #make sure this isn't too small
    mxrot = (0.5,)
    minwidth_rot = fill(0.002, 3)
    function worker(frame_n)
        @printf("===%d proc processing %d frame=====",  myid(), frame_n)
        @time tform, mm =qd_rigid(centered(img[:,:,1]), centered(img[:,:,frame_n]), 
        mxshift, mxrot, minwidth_rot, SD; maxevals=1000, rtol=0, fvalue=0.0002)
        return tform, mm
    end
end

time_step = size(img,3)
@time result = pmap(worker, 2:time_step; retry_delays = zeros(3), on_error=ex->(AffineMap(fill(NaN, (2,2)), fill(NaN, 2)),NaN)) #time_step

using JLD2, FileIO
array_result = Any[]
for (i, item) in enumerate(result)
    try
        push!(array_result, (item[1].linear, item[1].translation, item[2]))
        if isnan(item[2])
            print("Found NaN result in ", i)
        end
    catch
        print("Found NaN result in ", i)
        push!(array_result, NaN)
    end
end
@save "parallel_reg_result_array.jld2" array_result
using Serialization
serialize(open("parallel_reg_result.jls","w"),result) # Use IO stream instead of others
# convert the dataformat to more savable one!
