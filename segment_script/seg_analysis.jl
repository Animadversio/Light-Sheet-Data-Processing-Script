using CoordinateTransformations
using CellSegmentation, Images, FakeCells, Test, FileIO, JLD2, Statistics
using TileTrees, BoxTrees, IntervalSets
using LinearAlgebra
using Serialization

ttree = deserialize("/Users/binxu/Holy_Optical_Imaging/segment_script/segment_ttree.jls")
using Plots
pyplot()
plot(ttree.tiles[1].T)
