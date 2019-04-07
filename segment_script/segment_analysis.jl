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
