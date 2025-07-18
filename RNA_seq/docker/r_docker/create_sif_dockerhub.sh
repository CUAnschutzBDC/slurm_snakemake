#!/bin/bash 

#SBATCH --job-name=sif
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem=6gb
#SBATCH --output=logs/sif.out
#SBATCH --partition=amilan

singularity pull --name rnaseq_rv2_dockerhub.sif docker://kwellswrasman/rna_seq_r:v2