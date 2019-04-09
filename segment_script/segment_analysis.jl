using TileTreesGUI
using Serialization, FileIO, JLD2
using BoxTrees, TileTrees, TileTreesGUI, ImageView
using Images
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
f = load(tree_path)
ttree, reftree = f["ttree"], f["reftree"];

ttree = deserialize("/Volumes/Seagate_Backup_Binxu/stitched1/segment_ttreenew.jls")

img = load(img_path);
tileaction = initialize_triage_action(ttree);
triage_gui(ttree, tileaction, img)

## load data
img = load(img_path);
tree_path = "D:\\stitched2\\segment_ttree.jls"
ttree2 = deserialize(tree_path)
tileaction2 = initialize_triage_action(ttree2);
triage_gui(ttree2, tileaction2, img)
# %%
# using CellSegmentation, MergePairwise, TiledFactorizations
# ttreenew = triage_actions!(deepcopy(ttree), tileaction)
# tiled_nf!(ttreenew, reftree, initT=false)
# output_dir = "D:\\stitched2\\"
# serialize(open(joinpath(output_dir,"segment_ttreenew.jls"),"w"),ttreenew)
#
# triage_gui(ttreenew, tileaction, img, reftree)

#---
using Plots
#---
spatial_sum = [sum(tile.S) for tile in ttree2.tiles];
spatial_max = [maximum(tile.S) for tile in ttree2.tiles];
temporal_sum = [sum(tile.T) for tile in ttree2.tiles];
temporal_max = [maximum(tile.T) for tile in ttree2.tiles];

histogram(map(log, temporal_sum), )

scatter(spatial_sum, spatial_max,s=1000*temporal_max);
xlabel!("spatial_sum");
ylabel!("spatial_max");

scatter(map(log, temporal_sum),spatial_sum);xlabel!("log(temporal_sum)");ylabel!("spatial_sum")
scatter(map(log, temporal_max),spatial_sum);xlabel!("log(temporal_max)");ylabel!("spatial_sum")
scatter(map(log, temporal_max),spatial_sum, zcolor=spatial_max, alpha=0.5);xlabel!("log(temporal_max)");ylabel!("spatial_max")

using FFTW, Statistics
tile_id = 1051
T_power = map(abs2,FFTW.fft(ttree2.tiles[tile_id].T .- mean(ttree2.tiles[tile_id].T)));
plot(T_power)
ylims!(-0.05,1)
# Good 1005 1010 944
# Bad 1007
#%% Population distribution analysis
@time begin
temporal_power_arr = zeros(Float32,(length(ttree2.tiles), length(ttree2.tiles[1].T)))
for tile_id=1:length(ttree2.tiles)
    T_power = map(abs2,FFTW.fft(ttree2.tiles[tile_id].T .- mean(ttree2.tiles[tile_id].T)));
    temporal_power_arr[tile_id, :] = T_power
end
end
L0band = sum(temporal_power_arr[:,50:500], dims=2);
L1band = sum(temporal_power_arr[:,500:2000], dims=2);
Mband = sum(temporal_power_arr[:,2000:5000], dims=2);
Hband = sum(temporal_power_arr[:,5000:10000], dims=2);

scatter(log.(L0band), log.(L1band), xlabel="L0band",ylabel="L1band")
scatter(log.(L0band), log.(Mband), xlabel="L0band",ylabel="Mband")
scatter(log.(L0band), log.(Mband), xlabel="L0band",ylabel="Mband")

# --- Build mask from features values
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

# ---- Low Pass Filtering of traces
using DSP
trace_id = 946
responsetype = Lowpass(0.5; fs=20);
designmethod = Butterworth(4);
trace_filt = filt(digitalfilter(responsetype, designmethod), ttree2.tiles[trace_id].T);
plot(ttree2.tiles[trace_id].T)
plot!(trace_filt)
