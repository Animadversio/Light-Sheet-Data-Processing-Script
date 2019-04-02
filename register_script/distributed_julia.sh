#PBS -N Distrib_registration_100hz_2min_fish2
# Pretest of the ability of new model on new aligned volume

#PBS -l nodes=25:ppn=8,walltime=4:00:00,mem=30gb
# Request less resource
# Specify the default queue for the fastest nodes
#PBS -m be
#PBS -q old

cat $PBS_NODEFILE
export nodes_name="`cat $PBS_NODEFILE`"
# "$(echo "`cat $PBS_NODEFILE`" | head -n 1 | tail -1)"
module load julia
julia ~/Holy_registeration/cluster_script.jl