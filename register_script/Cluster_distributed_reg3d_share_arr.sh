#PBS -N Distrib_registration3d_array
# Pretest of the ability of new model on new aligned volume
#PBS -l nodes=40:ppn=8,walltime=24:00:00,mem=40gb
# Request less resource
# Specify the default queue for the fastest nodes
#PBS -m be
#PBS -q old
#PBS -t 1
# cat $PBS_NODEFILE
export nodes_name="`cat $PBS_NODEFILE`"
# "$(echo "`cat $PBS_NODEFILE`" | head -n 1 | tail -1)"
img_names='stitched.nhdr'
output_names='stitched/'

"$(echo "$param_list" | head -n $PBS_ARRAYID | tail -1)"
export img_path="/scratch/binxu.wang/OCPI_data/"
img_path+="$(echo "$img_names" | head -n $PBS_ARRAYID | tail -1)"
export save_dir="/scratch/binxu.wang/OCPI_data/"
save_dir+="$(echo "$output_names" | head -n $PBS_ARRAYID | tail -1)"

echo $img_path
echo $save_dir

module load julia
julia ~/Holy_registeration/cluster_register3d_script_share_arr.jl
