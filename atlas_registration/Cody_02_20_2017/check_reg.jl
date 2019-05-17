using Images, ImageView#ImagePlayer
cd("/Users/binxu/Holy_Optical_Imaging/Registrate_Atlas/Cody_02_20_2017/")
img0name = "resampled/forebrain_ocpi_red.nrrd"
imgname = "forebrain_ocpi_red_fixed_refbrain_prealigned_moving_Warped.nrrd"

function check_reg(img0name, imgname)

    img = load(imgname, mmap=false)
    img0 = load(img0name, mmap=false)

    img = AxisArray(img.data./maximum(img), axes(img))
    img0 = AxisArray(img0.data./maximum(img0), axes(img0))

    imshow(colorview(RGB, img0.data, zeroarray, img.data))

    return img, img0
end


check_reg(img0name, imgname);
