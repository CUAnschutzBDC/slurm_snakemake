#!/bin/bash 

#SBATCH --job-name=snakemake
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=1gb
#SBATCH --output=logs/snakemake_%J.out
#SBATCH --partition=amilan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kristen.wells-wrasman@cuanschutz.edu

set -o nounset -o pipefail -o errexit -x

mkdir -p logs

module load singularity/3.7.4

profile="slurm_profiles/slurm_profile_singularity"

# Function to run snakemake
run_snakemake() {
    local num_jobs=$1
    local config_file=$2
    local profile=$3

    snakemake \
        --snakefile Snakefile \
        --jobs $num_jobs \
        --latency-wait 60 \
        --rerun-incomplete \
        --configfile $config_file \
        --profile $profile \
        --use-singularity \
        --restart-times 1
}

run_snakemake 12 config.yaml $profile