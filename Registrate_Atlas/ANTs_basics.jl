using Images, ImageView
using Images, Unitful, Interpolations, Rotations, CoordinateTransformations
using ANTsRegistration
if Sys.isapple()
    cd("/Users/binxu/Holy_Optical_Imaging/Registrate_Atlas/Cody_02_20_2017/")
elseif Sys.iswindows()
    cd(raw"D:\cross_modal_registration\huc_gc6f_vglut_lox_R_lox_G\both_preprocessed\02_20_2017")
end
refbrain_path = "resampled/refbrain.nrrd"
refbrain = load(refbrain_path);

using MAT
if Sys.isapple()
    cd("/Users/binxu/Holy_Optical_Imaging/matlab_visualization")
elseif Sys.iswindows()
    cd(raw"D:\Light-Sheet-Data-Processing-Script\matlab_visualization")
end
vol = matread("image_vol.mat")
image_vol = vol["image_vol"];
vol = nothing

# Translation([0,-130,0])

# https://julialang.org/blog/2017/04/offset-arrays
# Resize / resampling

# Original OCPI imaging spacing
ps = [0.65;0.65;5.0] * Unitful.μm # [pixelspacing(img0)...]
imgboxsz = ps.*([size(image_vol)...] .- 1)
new_spacing = [1.3;1.3;3.0] * Unitful.μm
newsz = map((x,y)->floor(Int,x/y), imgboxsz, new_spacing) .+ 1
# Spatial filtering
if any(map((x,y)->x>y, size(image_vol), newsz)) #see if we need to smooth to avoid aliasing
    print("Smoothing image to avoid aliasing...\n")
    sigmas = map((o,n)->0.75f0*o/n, size(image_vol), newsz)
    @show kern = KernelFactors.gaussian(sigmas)   # from ImageFiltering
    #imgr = imresize(imfilter(img, kern, NA()), sz)
    @inbounds filtered = imfilter(image_vol, kern, NA())
    #filtered = imfilter_gaussian(img0.data, sigmas)
else
    filtered = image_vol
end
imginterp = Interpolations.interpolate(filtered, BSpline(Linear()), OnGrid())
# get fraction of the original pixelspacing unit that we must offset before we start sampling
edge_lengths = imgboxsz .- (newsz .- 1).*new_spacing
sampling_start = (edge_lengths/2)./ps
interpidxs = map((x,y,z)->range(1.0+x, stop = y-x, length = z), sampling_start, size(image_vol), newsz)
newimg = imginterp[interpidxs...]
# show the moving volume after changing pixel size.
imshow(newimg)

##########################
# Manually initial registration and visualize
function replace_nan!(x::AbstractArray{T}) where T
    for i = eachindex(x)
        if isnan(x[i])
            x[i] = zero(T)
        end
    end
end
rot = RotZ(1.0*pi/2)#Rotation([0,0,1],-1.0*pi/2)
trans = Translation([330,-140,-29])
newimg_warp = warp(newimg, trans ∘ LinearMap(rot) ); # transform
replace_nan!(newimg_warp)
# Color Merging View
img_pv = imshow(colorview(RGB, paddedviews(0, refbrain./maximum(refbrain), zeroarray, newimg_warp./maximum(newimg_warp))...))
# 136, 168; 137 317
norm_refbrain = refbrain./maximum(refbrain);
norm_warp_ocpi = newimg_warp./maximum(newimg_warp);
# 237 136 to 203 172

# ANTsRegistration
new_ps = (1.3, 1.3, 3.0) .* Unitful.μm
using OffsetArrays
offset_ref = OffsetArray(norm_refbrain[120:340,120:340,20:100], 120:340,120:340,20:100);
fixed = AxisArray(norm_refbrain[120:340,120:340,20:100], (:x, :y, :z), new_ps);#norm_refbrain[120:340,120:340,20:100]
# fixed = AxisArray(norm_refbrain[120:340,120:340,20:100],   (:x, :y, :z), new_ps);
moving = AxisArray(norm_warp_ocpi, (:x, :y, :z), new_ps);
orig_mov = AxisArray(image_vol, (:x, :y, :z), (0.65, 0.65, 5.0) .* Unitful.μm);
# ps = (0.65, 0.65, 5.0) .* Unitful.μm
stage = Stage(fixed, Global("Rigid"), MI(), (4,2,1), (4u"μm",2u"μm",1u"μm"), (1000,500,100))
# Both of the following should register well: fixed0 and moving0 are isotropic, so
# the rotation is a rotation; fixed and moving encode their pixel spacing, so even
# though "squashed" horizontally they also differ by a rotation (once the pixel spacing is accounted for)
# rigid = Stage(fixed, Global("Rigid"))
syn = Stage(fixed, SyN())
img_warp0 = register(fixed, moving, [stage, syn];histmatch=true, verbose=true); # /255 , syn
# img_warp0 = register(fixed, orig_mov, [stage, syn];histmatch=true); # /255 , syn
img_ANTs_input = imshow(colorview(RGB, paddedviews(0, fixed, zeroarray, moving)...))
img_ANTs_syn = imshow(colorview(RGB, paddedviews(0, norm_refbrain[120:340,120:340,20:100], zeroarray, img_warp0)...))

###############################################################
#%% Register the reference brain together
using FixedPointNumbers
if Sys.isapple()
    cd("/Users/binxu/Holy_Optical_Imaging/Registrate_Atlas")
elseif Sys.iswindows()
    cd(raw"D:\Light-Sheet-Data-Processing-Script\Registrate_Atlas")
end
ref_HucGcamp = load("MPIN-Atlas__Reference_brains__Live__HuCGCaMP5G.nrrd");
ref_Huc = load("MPIN-Atlas__Reference_brains__Fixed__HuCnlsGCaMP.nrrd");
begin # print the volume info
    rng = map(ax->axisvalues(ax)[1],ref_Huc.axes);
    ps_Huc = map(step,rng);
    print("Huc fixed brain volume size\n")
    display(rng)
    print("Pixel size: ")
    display(ps_Huc)

    rng = map(ax->axisvalues(ax)[1],ref_HucGcamp.axes);
    ps_HucGcamp = map(step,rng);
    print("HucGCamp5 life brain volume size\n")
    display(rng)
    print("Pixel size: ")
    display(ps_HucGcamp)
end
# imshow(ref_Huc)
img_ref_match = imshow(colorview(RGB, paddedviews(0, reinterpret.(N0f8,ref_HucGcamp),
    zeroarray, reinterpret.(N0f8,ref_Huc))...))

# %% Resample
ps = pixelspacing(ref_HucGcamp) # [pixelspacing(img0)...]
imgboxsz = ps_HucGcamp.*([size(ref_HucGcamp)...] .- 1)
new_spacing = ps_Huc
newsz = map((x,y)->floor(Int,x/y), imgboxsz, new_spacing) .+ 1
# Spatial filtering
if any(map((x,y)->x>y, size(ref_HucGcamp), newsz)) #see if we need to smooth to avoid aliasing
    print("Smoothing image to avoid aliasing...\n")
    sigmas = map((o,n)->0.75f0*o/n, size(ref_HucGcamp), newsz)
    @show kern = KernelFactors.gaussian(sigmas)   # from ImageFiltering
    #imgr = imresize(imfilter(img, kern, NA()), sz)
    @inbounds filtered = imfilter(ref_HucGcamp, kern, NA())
    #filtered = imfilter_gaussian(img0.data, sigmas)
else
    filtered = image_vol
end
imginterp = Interpolations.interpolate(filtered, BSpline(Linear()), OnGrid())
# get fraction of the original pixelspacing unit that we must offset before we start sampling
edge_lengths = imgboxsz .- (newsz .- 1).*new_spacing;
sampling_start = (edge_lengths/2)./ps;
interpidxs = map((x,y,z)->range(1.0+x, stop = y-x, length = z), sampling_start, size(ref_HucGcamp), newsz);
ref_HucGcamp_resamp = imginterp[interpidxs...];
ref_HucGcamp_resamp = AxisArray(ref_HucGcamp_resamp, (:x, :y, :z), new_spacing) # add unit to axis
# %%
img_ref_match = imshow(colorview(RGB, paddedviews(0, reinterpret.(N0f8,ref_Huc),
    zeroarray, ref_HucGcamp_resamp)...))
# ANTsregister the life brain and the fixed brain
conv_data = convert(Array{Float32,3},fixed_ref.data);
fixed_ref = AxisArray(conv_data, (:x, :y, :z), new_spacing);#AxisArray(ref_Huc, (:x, :y, :z), new_ps);
# moving_ref = ref_HucGcamp_resamp;# ref_HucGcamp # AxisArray(ref_HucGcamp, (:x, :y, :z), new_ps);
moving_ref = AxisArray(convert(Array{Float32,3},ref_HucGcamp_resamp.data), (:x, :y, :z), new_spacing);# ref_HucGcamp # AxisArray(ref_HucGcamp, (:x, :y, :z), new_ps);
stage = Stage(fixed_ref, Global("Rigid"), MI(), (16,4,2,1), (16u"μm",4u"μm",2u"μm",1u"μm"), (2000,1000,500,100))
# syn = Stage(fixed, SyN())
ref_GCamp_warp = register(fixed_ref, moving_ref, [stage]; histmatch=true, verbose=true); # /255 , syn
# img_ANTs_input = imshow(colorview(RGB, paddedviews(0, fixed, zeroarray, moving)...))
