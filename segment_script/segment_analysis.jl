using TileTreesGUI
using Serialization, FileIO, JLD2
using BoxTrees, TileTrees, TileTreesGUI, ImageView
using Images
tree_path = "D:\\stitched1\\segment_tree.jld2"
f = load(tree_path)
ttree, reftree = f["ttree"], f["reftree"]

img_path = "F:\\stitched.nhdr"
img = load(img_path)
tileaction = initialize_triage_action(ttree);
triage_gui(ttree, tileaction, img)
