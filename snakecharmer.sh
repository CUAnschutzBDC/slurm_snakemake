#!/bin/bash 

#SBATCH --job-name=snakemake
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --mem=1gb
#SBATCH --output=logs/snakemake_%J.out
#SBATCH --partition=amilan


# Run snakemake pipeline
snakemake \
    --snakefile Snakefile \
    --jobs 15 \
    --latency-wait 60 \
    --ignore-incomplete \
    --profile profile/ 
