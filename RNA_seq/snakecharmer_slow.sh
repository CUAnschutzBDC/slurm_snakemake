#!/bin/bash 

#SBATCH --job-name=snakemake
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=60gb
#SBATCH --output=logs/snakemake_%J.out
#SBATCH --partition=amilan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kristen.wells-wrasman@cuanschutz.edu

set -o nounset -o pipefail -o errexit -x

mkdir -p logs
module load anaconda
conda activate snakemake8

#export APPTAINER_BIND="/scratch/alpine:/scratch/alpine,/pl/active/Anschutz_BDC:/pl/active/Anschutz_BDC"

# Run snakemake pipeline
snakemake \
    --snakefile Snakefile \
    --configfile config.yaml \
    --jobs 1 \
    --latency-wait 60 \
    --rerun-incomplete \
    --workflow-profile profiles/default \
    --executor local 