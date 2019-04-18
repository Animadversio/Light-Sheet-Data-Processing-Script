# Fish Visualization

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
#---
img = load(img_path);
using Makie
scene = Scene()
#plot the space inside
rx = 1:size(img,1) * 0.65
ry = 1:size(img,2) * 0.65
rz = 1:size(img,3) * 5
mat = 100*img.data[:,:,:,10000]#[(x.^2 + y.^2 + z.^2) for x = r, y = r, z = r]
mat2 = mat .* (mat .> 0.25)
volume(rx, ry, rz, 5* mat2, algorithm = :absorptionrgba)

scene2 = scene()
plotly()
histogram(collect(Iterators.flatten(mat)))
