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


tree_path = "D:\\stitched2\\segment_ttree.jls"
ttree2 = deserialize(tree_path)
tileaction2 = initialize_triage_action(ttree2);
triage_gui(ttree2, tileaction2, img)

using CellSegmentation, MergePairwise, TiledFactorizations
ttreenew = triage_actions!(deepcopy(ttree), tileaction)
tiled_nf!(ttreenew, reftree, initT=false)
output_dir = "D:\\stitched1\\"
serialize(open(joinpath(output_dir,"segment_ttreenew.jls"),"w"),ttreenew)

triage_gui(ttreenew, tileaction, img, reftree)

segment_prepare()

spatial_sum = [sum(tile.S) for tile in ttree2.tiles];
spatial_max = [maximum(tile.S) for tile in ttree2.tiles];
temporal_sum = [sum(tile.T) for tile in ttree2.tiles];
temporal_max = [maximum(tile.T) for tile in ttree2.tiles];
histogram(log(temporal_sum), )

scatter(spatial_sum, spatial_max,s=1000*temporal_max);
xlabel!("spatial_sum");
ylabel!("spatial_max");

scatter(map(log, temporal_sum),spatial_sum);xlabel!("log(temporal_sum)");ylabel!("spatial_max")
scatter(map(log, temporal_max),spatial_sum);xlabel!("log(temporal_max)");ylabel!("spatial_max")
scatter(map(log, temporal_max),spatial_sum, zcolor=spatial_max, alpha=0.5);xlabel!("log(temporal_max)");ylabel!("spatial_max")
