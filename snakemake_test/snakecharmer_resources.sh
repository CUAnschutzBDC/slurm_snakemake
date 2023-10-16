#!/bin/bash 

#SBATCH --job-name=snakemake
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --mem=1gb
#SBATCH --output=logs/snakemake_%J.out
#SBATCH --partition=amilan

profile="slurm_profiles/slurm_profile_resources" # Or set to slurm_profile_singularity if using singularity

mkdir -p logs

# Run snakemake pipeline
snakemake \
    --snakefile Snakefile_resources \
    --jobs 15 \
    --latency-wait 60 \
    --ignore-incomplete \
    --profile $profile
