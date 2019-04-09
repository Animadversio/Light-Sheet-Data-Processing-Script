#---
using TileTreesGUI
using Serialization, FileIO, JLD2
using BoxTrees, TileTrees, TileTreesGUI, ImageView
using Images

#---
if Sys.isapple()
    tree_path = "/Volumes/Seagate_Backup_Binxu/stitched1/segment_tree.jld2"
    img_path = "/Volumes/Seagate_Backup_Binxu/stitched.nhdr"
elseif Sys.iswindows()
    img_path = "D:\\stitched1\\stitched.nhdr"
    tree_path = "D:\\stitched1\\segment_tree.jld2"
elseif Sys.islinux()
    img_path = "/mnt/binxu001/stitched.nhdr"#"F:\\stitched.nhdr"
    tree_path = "/mnt/binxu001/stitched1/segment_tree.jld2"
end
#---
img = load(img_path);
tree_path = "D:\\stitched2\\segment_ttree.jls"
ttree2 = deserialize(tree_path);
tileaction2 = initialize_triage_action(ttree2);

#---
using CellSegmentation, MergePairwise, TiledFactorizations
using FFTW, Statistics
#--- Compute features
spatial_sum = [sum(tile.S) for tile in ttree2.tiles];
spatial_max = [maximum(tile.S) for tile in ttree2.tiles];
temporal_sum = [sum(tile.T) for tile in ttree2.tiles];
temporal_max = [maximum(tile.T) for tile in ttree2.tiles];
@time begin
temporal_power_arr = zeros(Float32,(length(ttree2.tiles), length(ttree2.tiles[1].T)))
for tile_id=1:length(ttree2.tiles)
    trace = ttree2.tiles[tile_id].T;
    T_power = abs2.(FFTW.fft((trace .- mean(trace))./std(trace)))/length(trace);
    temporal_power_arr[tile_id, :] = T_power
end
end
Lband = sum(temporal_power_arr[:,10:20], dims=2);
L0band = sum(temporal_power_arr[:,20:50], dims=2);
L1band = sum(temporal_power_arr[:,50:100], dims=2);
L2band = sum(temporal_power_arr[:,100:200], dims=2);
L3band = sum(temporal_power_arr[:,200:500], dims=2);
L4band = sum(temporal_power_arr[:,500:1000], dims=2);
L5band = sum(temporal_power_arr[:,1000:2000], dims=2);
Mband = sum(temporal_power_arr[:,2000:5000], dims=2);
Hband = sum(temporal_power_arr[:,5000:10000], dims=2);
#---
using Plots #,Plotly,GR
# plotly()
# GR()
pyplot()
scatter(log.(temporal_sum),spatial_sum);xlabel!("log(temporal_sum)");ylabel!("spatial_sum")
scatter(log.(temporal_max),spatial_sum);xlabel!("log(temporal_max)");ylabel!("spatial_sum")
scatter(log.(temporal_max),spatial_sum, zcolor=spatial_max, alpha=0.5);xlabel!("log(temporal_max)");ylabel!("spatial_max")

scatter(log.(Lband), log.(Mband), xlabel="Lband",ylabel="Mband")
scatter(log.(Lband), log.(L1band), xlabel="Lband",ylabel="L1band")
scatter((Lband), (L1band), xlabel="Lband",ylabel="L1band")
scatter(log.(L0band), log.(L1band), xlabel="L0band",ylabel="L1band")
scatter(log.(L0band), log.(Mband), xlabel="L0band",ylabel="Mband")
scatter(log.(L0band), log.(Mband), xlabel="L0band",ylabel="Mband")

#%%

function find_nonzero(array)
    indx = [];
    for i = 1:length(array)
        if array[i]!=0
        push!(indx,i)
        end
    end
return indx
end
tile_mask = (L0band.>3) .& (spatial_sum .> 10) .& (map(log, temporal_max) .> -6.5);
tile_idx  = find_nonzero(tile_mask);
# ---
function mask_action!(tile_mask, tileaction)
    for i = 1:length(tile_mask)
        if tile_mask[i]
            tileaction[i] = ("keep",)
        else
            tileaction[i] = ("delete",)
        end
    end
end
mask_action!(tile_mask, tileaction2)
ttree2_new = triage_actions!(deepcopy(ttree2), tileaction2)
# tiled_nf!(ttreenew, reftree, initT=false)
triage_gui(ttree2_new, tileaction2, img)


#%%
# ---- Low Pass Filtering of traces
using DSP
trace_id = 946
responsetype = Lowpass(0.5; fs=20);
designmethod = Butterworth(4);
trace_filt = filt(digitalfilter(responsetype, designmethod), ttree2.tiles[trace_id].T);
plot(ttree2.tiles[trace_id].T)
plot!(trace_filt)

using Printf
# --- View features of a tile
sampling_rate = 20;
tile_id = 1005
trace = ttree2.tiles[tile_id].T;
L = length(trace)
f_space = sampling_rate/2 * LinRange(0,1,ceil(Int,L/2+1));
T_power = abs2.(FFTW.fft((trace .- mean(trace))./std(trace)))/length(trace);
@printf("===Tile %d===\n", tile_id)
@printf("Sum spatial val: %.3f \t Max spatial val: %.3f
        Sum temporal val: %.3f \t Max temporal val: %.3f
        Centralized Band Power: [10,20]: %.3f \t [20,50]: %.3f \t
        [50,100]: %.3f \t  [100,200]: %.3f \t  [200,500]: %.3f \t  [500,1000]: %.3f \t
        [1000,2000]: %.3f \t [2000,5000]: %.3f \t [5000,10000]: %.3f \n",
        spatial_sum[tile_id], spatial_max[tile_id],
        temporal_sum[tile_id], temporal_max[tile_id], Lband[tile_id], L0band[tile_id],
        L1band[tile_id], L2band[tile_id], L3band[tile_id], L4band[tile_id],
        L5band[tile_id], Mband[tile_id], Hband[tile_id])

plot(f_space, T_power[1:ceil(Int,L/2+1)], label=string(tile_id), title="Spectrum Power of Tile",
    xlabel="freq(Hz)")
ylims!(-min(std(T_power),0.05), 3*std(T_power))
# 946 947 954 957 958 is really good tile! 955 weak cell but less SNR
# 981 is an distinct example
# 1003 is definitely bad
# 1004 is ful of noise
