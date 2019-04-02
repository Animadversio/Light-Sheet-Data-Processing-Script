#PBS -N Apply_Reg_Cell_Segmentation_Array
# Pretest of the ability of new model on new aligned volume
#PBS -l nodes=1:ppn=8,walltime=24:00:00,mem=45gb
# Request less resource
# Specify the default queue for the fastest nodes
#PBS -m be
#PBS -q old
#PBS -t 1,3,5

# cat $PBS_NODEFILE
# export nodes_name="`cat $PBS_NODEFILE`"
# "$(echo "`cat $PBS_NODEFILE`" | head -n 1 | tail -1)"
export nodes_name="`cat $PBS_NODEFILE`"
# "$(echo "`cat $PBS_NODEFILE`" | head -n 1 | tail -1)"
img_names='single_plane_100hz_5min.imagine
single_plane_100hz_2min_fish2_after_tricaine.imagine
single_plane_100hz_2min_fish2_after_tricaine2.imagine
single_plane_100hz_2min_fish2_after_tricaine3.imagine
single_plane_100hz_2min_fish2_after_tricaine4.imagine'
output_names='single_plane_100hz_5min/
single_plane_100hz_2min_fish2_after_tricaine/
single_plane_100hz_2min_fish2_after_tricaine2/
single_plane_100hz_2min_fish2_after_tricaine3/
single_plane_100hz_2min_fish2_after_tricaine4/'
# "$(echo "$param_list" | head -n $PBS_ARRAYID | tail -1)"
export img_path="/scratch/binxu.wang/OCPI_data/"
img_path+="$(echo "$img_names" | head -n $PBS_ARRAYID | tail -1)"
export save_dir="/scratch/binxu.wang/OCPI_data/"
save_dir+="$(echo "$output_names" | head -n $PBS_ARRAYID | tail -1)"

echo $img_path
echo $save_dir
module load julia
julia ~/Holy_registeration/apply_registration.jl
julia ~/Holy_segmentation/cell_segmentation_script.jl
