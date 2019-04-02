#PBS -N Distrib_registration_100hz_array
# Pretest of the ability of new model on new aligned volume
#PBS -l nodes=30:ppn=8,walltime=12:00:00,mem=40gb
# Request less resource
# Specify the default queue for the fastest nodes
#PBS -m be
#PBS -q old
#PBS -t 7
# cat $PBS_NODEFILE
export nodes_name="`cat $PBS_NODEFILE`"
# "$(echo "`cat $PBS_NODEFILE`" | head -n 1 | tail -1)"
img_names='single_plane_100hz_5min.imagine
single_plane_100hz_2min_fish2_after_tricaine.imagine
single_plane_100hz_2min_fish2_after_tricaine2.imagine
single_plane_100hz_2min_fish2_after_tricaine3.imagine
single_plane_100hz_2min_fish2_after_tricaine4.imagine
03_26_2019/40x_4dpf_fish1_gcamp6f.imagine
03_27_2019/dpf5_fish2_gfp_2d_dynamic.imagine'

output_names='single_plane_100hz_5min/
single_plane_100hz_2min_fish2_after_tricaine/
single_plane_100hz_2min_fish2_after_tricaine2/
single_plane_100hz_2min_fish2_after_tricaine3/
single_plane_100hz_2min_fish2_after_tricaine4/
03_26_2019/40x_4dpf_fish1_gcamp6f/
03_27_2019/dpf5_fish2_gfp_2d_dynamic/'

# "$(echo "$param_list" | head -n $PBS_ARRAYID | tail -1)"
export img_path="/scratch/binxu.wang/OCPI_data/"
img_path+="$(echo "$img_names" | head -n $PBS_ARRAYID | tail -1)"
export save_dir="/scratch/binxu.wang/OCPI_data/"
save_dir+="$(echo "$output_names" | head -n $PBS_ARRAYID | tail -1)"

echo $img_path
echo $save_dir

module load julia
julia ~/Holy_registeration/cluster_register_script_share_arr.jl
