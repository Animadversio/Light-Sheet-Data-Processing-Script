using Images, ImageView
cd("/Users/binxu/Holy_Optical_Imaging/Registrate_Atlas/Cody_02_20_2017/")
refbrain_path = "resampled/refbrain.nrrd"
refbrain = load(refbrain_path)
imshow(refbrain)
using MAT
cd("/Users/binxu/Holy_Optical_Imaging/matlab_visualization")
vol = matread("image_vol.mat")
image_vol = vol["image_vol"];
vol = nothing

#%%
