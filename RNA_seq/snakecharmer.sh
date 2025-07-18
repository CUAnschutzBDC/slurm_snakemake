#!/bin/bash 

#SBATCH --job-name=snakemake
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=1gb
#SBATCH --output=logs/snakemake_%J.out
#SBATCH --partition=amilan
#SBATCH --mail-type=ALL
#SBATCH --qos=normal
#SBATCH --mail-user=kristen.wells-wrasman@cuanschutz.edu

set -o nounset -o pipefail -o errexit -x

mkdir -p logs
module load anaconda
conda activate snakemake8

# For some reason this is set in acompile and vscode but not on login nodes
# This line fixes the error 
# srun: fatal: SLURM_MEM_PER_CPU, SLURM_MEM_PER_GPU, and SLURM_MEM_PER_NODE are mutually exclusive.
# in spaned snakemake jobs 
unset SLURM_MEM_PER_CPU

#export APPTAINER_BIND="/scratch/alpine:/scratch/alpine,/pl/active/Anschutz_BDC:/pl/active/Anschutz_BDC"

# Run snakemake pipeline
snakemake \
    --snakefile Snakefile \
    --configfile config.yaml \
    --jobs 12 \
    --latency-wait 60 \
    --rerun-incomplete \
    --workflow-profile profiles/default 