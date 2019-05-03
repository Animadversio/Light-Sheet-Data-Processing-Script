using NIfTI, Images, ImageView
using ImageCore, Colors

ni = niread("D:\\cross_modal_registration\\huc_gc6f_vglut_lox_R_lox_G\\both_preprocessed\\02_20_
2017\\forebrain_ocpi_red_fixed_refbrain_prealigned_moving1InverseWarp.nii.gz")# , mmap=true
# forebrain_ocpi_red_fixed_refbrain_prealigned_moving1InverseWarp.nii.gz
# a 5 dimension matrix, seems to be the warp vector field
# imshow(ni[:,:,:,1,:])
ni_warp = niread(raw"D:\cross_modal_registration\huc_gc6f_vglut_lox_R_lox_G\both_preprocessed\02_20_2017\forebrain_ocpi_red_fixed_refbrain_prealigned_moving1Warp.nii.gz");
# Visualize the warp field as a color tensor!
rgb_vol = colorview(RGB, ni_warp[:,:,:,1,1], ni_warp[:,:,:,1,2], ni_warp[:,:,:,1,3]);
rgb_vol_in = colorview(RGB, ni[:,:,:,1,1], ni[:,:,:,1,2], ni[:,:,:,1,3]);
imshow(rgb_vol_in)
imshow(rgb_vol)


#%%
img = load(raw"D:\cross_modal_registration\huc_gc6f_vglut_lox_R_lox_G\both_preprocessed\02_20_2017\forebrain_ocpi_red_fixed_refbrain_prealigned_moving_Warped.nrrd")
imshow(img)

using MAT
f= matopen(raw"D:\cross_modal_registration\huc_gc6f_vglut_lox_R_lox_G\both_preprocessed\02_20_2017\testdir\testdirforebrain_ocpi_red_fixed_refbrain_prealigned_moving0GenericAffine.mat")
# %%
img = load(raw"D:\cross_modal_registration\huc_gc6f_vglut_lox_R_lox_G\both_preprocessed\02_20_2017\resampled\forebrain_ocpi_green.nrrd")
imshow(img)

img = load(raw"D:\cross_modal_registration\huc_gc6f_vglut_lox_R_lox_G\both_preprocessed\02_20_2017\resampled\refbrain_prealigned.nrrd")
imshow(img)

img = load(raw"D:\cross_modal_registration\huc_gc6f_vglut_lox_R_lox_G\both_preprocessed\02_20_2017\resampled\whole_ocpi_red.nrrd")
imshow(img)

habenula_img = load(raw"D:\cross_modal_registration\huc_gc6f_vglut_lox_R_lox_G\both_preprocessed\02_20_2017\atlas_warped_for_this_fish\anatomy\habenula.nrrd")
imshow(habenula_img)

ref_img = load(raw"D:\cross_modal_registration\huc_gc6f_vglut_lox_R_lox_G\confocal\december_2016\pre_aligned\refbrain_prealigned_f32_no_nans.nrrd")
imshow(habenula_img)
