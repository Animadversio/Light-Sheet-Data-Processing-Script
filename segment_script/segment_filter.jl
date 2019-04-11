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
if Sys.isapple()
    tree_path = "/Users/binxu/Holy_Optical_Imaging/segment_ttree.jls"
elseif Sys.iswindows()
    tree_path = "D:\\stitched2\\segment_ttree.jls"
end
ttree2 = deserialize(tree_path);
tileaction2 = initialize_triage_action(ttree2);

#---
using CellSegmentation, MergePairwise, TiledFactorizations
using FFTW, Statistics, DSP, StatsBase
using DSP
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
LLband = sum(temporal_power_arr[:,1:10], dims=2);
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
using Plots, Printf #,Plotly,GR
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
# ---- Low Pass Filtering of traces

trace_id = 946
responsetype = Lowpass(0.5; fs=20);
designmethod = Butterworth(4);
trace_filt = filt(digitalfilter(responsetype, designmethod), ttree2.tiles[trace_id].T);
plot(ttree2.tiles[trace_id].T)
plot!(trace_filt)


# --- View features of a tile

f = Plots.font("DejaVu Sans", 12)
default(size=(1200,1000), guidefont=f, xtickfont=f, ytickfont=f, titlefont=font(14), legendfont=f)
function vis_calc_feature(tile_id, do_plot::Bool=true)
    sampling_rate = 20;
    trace = ttree2.tiles[tile_id].T;
    L = length(trace)
    f_space = sampling_rate/2 * LinRange(0,1,ceil(Int,L/2+1));
    T_power = abs2.(FFTW.fft((trace .- mean(trace))./std(trace)))/length(trace);

    feat_stats = [spatial_sum[tile_id], spatial_max[tile_id],
    temporal_sum[tile_id], temporal_max[tile_id],
    std(trace), skewness(trace), kurtosis(trace),
    Lband[tile_id], L0band[tile_id], L1band[tile_id], L2band[tile_id],
    L3band[tile_id], L4band[tile_id], L5band[tile_id], Mband[tile_id],
    Hband[tile_id]]
    @printf("===Tile %d===\n", tile_id)
    @printf("Sum spatial val: %.3f \t Max spatial val: %.3f
    Sum temporal val: %.3f \t Max temporal val: %.3f
    Trace std: %.3f \t Trace skewness: %.3f \t Trace kurtosis: %.3f \t
    Centralized Band Power: [10,20]: %.3f \t [20,50]: %.3f \t
    [50,100]: %.3f \t  [100,200]: %.3f \t  [200,500]: %.3f \t  [500,1000]: %.3f \t
    [1000,2000]: %.3f \t [2000,5000]: %.3f \t [5000,10000]: %.3f \n",
    feat_stats...)
    if !do_plot
        return feat_stats
    else
        plt0 = histogram(trace, bin=50, orientation=:horizontal, label=string(tile_id),
            alpha=0.8, legend=false, yticks=nothing, xticks=nothing)
        title!(@sprintf("Signal histogram"))
        Plots.annotate!(0.6*xlims()[2], 0.7*ylims()[2],
            text(@sprintf("std: %.3f \n skew: %.3f \n kurt: %.3f", feat_stats[5:7]...), 12))

        plt1 = plot(f_space, T_power[1:ceil(Int,L/2+1)], label=string(tile_id),
            xlabel="Freq(Hz)", legend=:bottom)
        xticks!(0:0.5:sampling_rate/2)
        vline!(collect(f_space)[[10,20,50,100,200,500,1000,2000,5000,10000]],
            label="Band seperation",alpha=0.8,linewidth=0.4, color=:red)
        title!(@sprintf("Spectrum Power of Tile %d", tile_id))
        ylims!(-min(std(T_power),0.05), 2 * percentile(Tpower, 99.5))#3*std(T_power))
        Plots.annotate!(sampling_rate/4, 0.9*ylims()[2],
            text(@sprintf("Band Power: [10,20]: %.3f [20,50]: %.3f
        [50,100]: %.3f [100,200]: %.3f  [200,500]: %.3f  [500,1000]: %.3f
        [1000,2000]: %.3f [2000,5000]: %.3f [5000,10000]: %.3f \n",feat_stats[8:end]...), 12))

        plt2 = plot((1:L)./sampling_rate, trace, xlabel="Time(s)", label=string(tile_id),
            linewidth=0.8, alpha=0.8)
        # plt_syn = plot(plt1, plt2, layout=grid(2,1,heights=[0.7,0.3]))
        htm1 = heatmap(ttree2.tiles[tile_id].S[:,:,1],aspect_ratio=:equal,grid=false,ticks=false,legend = :none)
        htm2 = heatmap(ttree2.tiles[tile_id].S[:,:,2],aspect_ratio=:equal,grid=false,ticks=false,legend = :none)
        htm3 = heatmap(ttree2.tiles[tile_id].S[:,:,3],aspect_ratio=:equal,grid=false,ticks=false,legend = :none)
        blk = plot(legend=false,grid=false,ticks=false,foreground_color_subplot=:white)
        l = @layout [  a{0.6h};
                    grid(1,2, widths=[0.8,0.2]);
                    grid(1,3)] # , widths=[0.2,0.2,0.2,0.2,0.2]
        plt_syn = plot(plt1, plt2, plt0, htm1, htm2, htm3,  layout=l)
        return feat_stats, plt_syn, [plt0,plt1,plt2]
    end
end
tile_id = 1320
trace = ttree2.tiles[tile_id].T;
stats, pltsyn, plt_arr = vis_calc_feature(tile_id)
display(pltsyn)
# 946 947 954 957 958 is really good tile! 955 weak cell but less SNR
# 981 is an distinct example
# 1003 is definitely bad
# 1004 1022 1023 is ful of noise.  1020 different kind of noise
# 1081 is purely noise
# 1025 is an example of long term ramping noise! Note the difference
@time for tile_id=tile_idx
    stats, pltsyn, plt_arr = vis_calc_feature(tile_id)
    savefig(pltsyn, "D:\\Holy Lab\\Proc_figures\\"*@sprintf("stats_%04d.png",tile_id))
end # 2303.141171 seconds (819.15 M allocations: 41.391 GiB, 0.82% gc time)
stats_arr = zeros(Float32,(length(ttree2.tiles), 16))
@time for tile_id=1:length(ttree2.tiles)
    stats = vis_calc_feature(tile_id, false);
    stats_arr[tile_id, :] = stats;
end # 42.071125 seconds (4.10 M allocations: 2.897 GiB, 1.12% gc time)
scatter(log.(stats_arr[:,9]),log.(abs.(stats_arr[:,7])))
# ---
function find_nonzero(array)
 indx = [];
 for i = 1:length(array)
     if array[i]!=0
     push!(indx,i)
     end
 end
return indx
end
tile_mask = (L0band.>3) .& (spatial_sum .> 10) .& (map(log, temporal_max) .> -6.5) .&
     (stats_arr[:,7] .> 2);
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


#%% Heatmap of population activity
@time begin
events_arr = zeros(Bool,(length(tile_idx), length(trace)))
for (idx, tile_id) in enumerate(tile_idx)
    trace = ttree2.tiles[tile_id].T
    events_arr[idx, :] = trace.>2*std(trace)
end
end #  1.225967 seconds (45.73 k allocations: 73.329 MiB, 0.53% gc time)
@time htmap = heatmap((1:length(trace))/20,1:length(tile_idx),events_arr,size=(1800,600),legend=:none,
    title="Events raster of >2 std activity", show=false);
savefig(htmap, "D:\\Holy Lab\\events_raster2.png")

#%% zscore map
@time begin
zscore_arr = zeros(Float32,(length(tile_idx), length(ttree2.tiles[1].T)))
for (idx, tile_id) in enumerate(tile_idx)
    trace = ttree2.tiles[tile_id].T
    zscore_arr[idx, :] = zscore(trace)
end
end #  1.225967 seconds (45.73 k allocations: 73.329 MiB, 0.53% gc time)

tmp = heatmap(zscore_arr[1:500,1:2000],size=(800,600),c=:OrRd,clims=(0,5),
    title="Z score activity");
savefig(tmp, "D:\\Holy Lab\\tmp.png")

@time htmap2 = heatmap((1:length(trace))/20,1:length(tile_idx),
    zscore_arr,c=:OrRd,clims=(0,5),size=(1800,600),
    title="Z score activity", show=false);
savefig(htmap2, "D:\\Holy Lab\\zscore_plot.png")
display(htmap2)

#%%
sumactivity = sum(zscore_arr[:,17000:20000], 1)
max_val, rel_idx = findmax(sumactivity)
synfire_idx = findall(mean(zscore_arr[:,18600:18620],dims=2).>5)
synfire_glob_idx = tile_idx[synfire_idx]
@time plt = plot((18200:19500)/20,
    transpose(zscore_arr[[coord.I[1] for coord in synfire_idx],18200:19500]),
    legend=:none,size=(1500,600),alpha=0.6,)
xlabel!("Time(s)")
ylabel!("z-score")
@time savefig(plt,"D:\\Holy Lab\\zscore_zoomin_plot.png")
# [coord.I[1] for coord in synfire_idx]

sumactivity = sum(zscore_arr[:,6200:8000], dims=1)
max_val, rel_idx = findmax(sumactivity)
synfire_idx = findall(mean(zscore_arr[:,7135:7145],dims=2).>5)
synfire_glob_idx = tile_idx[synfire_idx]
time_window = 7000:8000
@time plt = plot(time_window/20,
    transpose(zscore_arr[[coord.I[1] for coord in synfire_idx],time_window]),
    legend=:none,size=(1500,600),alpha=0.6,)
xlabel!("Time(s)")
ylabel!("z-score")
@time savefig(plt,"D:\\Holy Lab\\zscore_zoomin_plot2.png")



#%%% Save and output data!
using MAT
f = matopen("D:\\Holy Lab\\tile_traces.mat", "w")
write(f, "zscore_arr", zscore_arr)
write(f, "tile_idx", tile_idx)
write(f, "stats_arr", stats_arr)
close(f)


raw_traces = [ttree2.tiles[tile_id].T for tile_id in 1:length(ttree2.tiles)]
raw_traces2 = hcat(raw_traces...)
f = matopen("D:\\Holy Lab\\raw_traces.mat", "w")
write(f, "raw_traces", raw_traces2)
close(f)

tile_space = [ttree2.tiles[tile_id].S for tile_id in 1:length(ttree2.tiles)]
f = matopen("D:\\Holy Lab\\tile_space.mat", "w")
write(f, "tile_space", tile_space)
close(f)

tile_pos = [[[ttree2.tiles[tile_id].box.intervals[1].left, ttree2.tiles[tile_id].box.intervals[1].right],
             [ttree2.tiles[tile_id].box.intervals[2].left, ttree2.tiles[tile_id].box.intervals[2].right],
             [ttree2.tiles[tile_id].box.intervals[3].left, ttree2.tiles[tile_id].box.intervals[3].right]]
             for tile_id in 1:length(ttree2.tiles)]
f = matopen("D:\\Holy Lab\\tile_pos2.mat", "w")
write(f, "tile_pos", tile_pos)
close(f)

using NPZ
npzwrite("D:\\Holy Lab\\ttree_data.npz", Dict("zscore_arr" => zscore_arr, "tile_pos" => tile_pos, "raw_traces" => raw_traces,
                "tile_space" => tile_space, "tile_idx" => tile_idx))
