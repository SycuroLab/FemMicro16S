#!/bin/bash

#SBATCH --partition=synergy,cpu2019,cpu2021,cpu2022,cpu2023,cpu2017-bf05,cpu2019-bf05
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH --mem=38G
#SBATCH --error=record_dada2.%J.err
#SBATCH --output=record_dada2.%J.out

log_dir="$(pwd)"
log_file="logs/dada2-analysis.log.txt"
num_jobs=60

# Path to your Apptainer image with Snakemake installed
IMAGE="apptainer/snakemake-1.0.0.sif"

# Snakemake command wrapper (runs inside the .sif)
SNAKEMAKE="apptainer exec --bind $PWD:$PWD --pwd $PWD $IMAGE snakemake"

echo "Started at: $(date)"

# ---- Run Snakemake inside the container ----

# Unlock, in case a previous run was interrupted
$SNAKEMAKE --unlock

# Main Snakemake call
$SNAKEMAKE --singularity-args "-B /bulk" --rerun-triggers mtime --latency-wait 60 --rerun-incomplete  --cluster-config cluster.json --cluster 'sbatch --partition={cluster.partition} --cpus-per-task={cluster.cpus-per-task} --nodes={cluster.nodes} --ntasks={cluster.ntasks} --time={cluster.time} --mem={cluster.mem} --output={cluster.output} --error={cluster.error}' --jobs $num_jobs --use-conda --use-singularity &>> $log_dir/$log_file

output_dir=$(grep "output_dir" < config.yaml | cut -d ' ' -f2 | sed 's/"//g')
list_files=$(grep "sampletable" < config.yaml | cut -d ' ' -f2 | sed 's/"//g')

#Copying all snakemake/log files of the run in the output folder
snakemake_file_dir="${output_dir}/snakemake_files"
mkdir -p $snakemake_file_dir

cp $list_files $snakemake_file_dir
cp Snakefile $snakemake_file_dir
cp config.yaml $snakemake_file_dir
cp cluster.json $snakemake_file_dir
cp dada2_sbatch.sh $snakemake_file_dir 
cp -rf logs $snakemake_file_dir
cp -rf utils $snakemake_file_dir

bash Version_check.sh > used_tools_versions.txt

echo "Finished at: `date`"
