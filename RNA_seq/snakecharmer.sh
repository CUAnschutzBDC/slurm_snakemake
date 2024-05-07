#!/bin/bash 

#SBATCH --job-name=snakemake
#SBATCH --ntasks=1
#SBATCH --time=23:50:00
#SBATCH --mem=1gb
#SBATCH --output=logs/snakemake_%J.out
#SBATCH --partition=amilan

set -o nounset -o pipefail -o errexit -x

mkdir -p logs
module load anaconda
conda activate snakemake

#module load singularity/3.7.4

profile="slurm_profiles/slurm_profile_singularity"

# Run snakemake pipeline
snakemake \
    --snakefile Snakefile \
    --configfile config.yaml \
    --jobs 12 \
    --latency-wait 60 \
    --rerun-incomplete \
    --profile $profile \
    --use-singularity \
    --singularity-args "--bind /scratch/alpine/$USER --bind /pl/active/Anschutz_BDC --bind /tmp:/tmp"
    