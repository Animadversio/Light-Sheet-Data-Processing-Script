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


# --- View features of a tile

f = Plots.font("DejaVu Sans", 12)
default(size=(1200,600), guidefont=f, xtickfont=f, ytickfont=f, titlefont=font(14), legendfont=f)
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

        l = @layout [  a{0.7h};
                    grid(1,2, widths=[0.8,0.2])]
        plt_syn = plot(plt1, plt2, plt0, layout=l)
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
end
stats_arr = zeros()
@time for tile_id=1:length(ttree2.tiles)
    stats = vis_calc_feature(tile_id, false)
