using JLD2, FileIO
using Images
using CoordinateTransformations
using LinearAlgebra
using Serialization
result = deserialize(open("/Users/binxu/Holy_Optical_Imaging/register_script/parallel_reg_result.jls","r"))
# f = load("/Users/binxu/Holy_Optical_Imaging/register_script/parallel_reg_result_array.jld2")
AffineMap(reshape(collect(f["result"][1][1].linear.data),(2,2)),collect(f["result"][1][1].translation.data))

##
using Plots
# using PyPlot
plotly() # or pyplot()
plot([item[2] for item in result],linewidth=2,title="Minimum Loss")
plotly()
x_move = [item[1].translation[1] for item in result]
y_move = [item[1].translation[2] for item in result]
plot([x_move, y_move],linewidth=2,title="Translation vector",label=["x_mov", "y_mov"])
xlabel("Frame number")
legend()

angle = [item[1].linear[1,2] for item in result]
plot(angle,linewidth=2,title="Rotation Angle",label=["Angle (rad)"])
xlabel("Frame number")
