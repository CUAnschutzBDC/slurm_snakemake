#!/bin/bash 

#SBATCH --job-name=download_fastq
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem=1gb
#SBATCH --output=logs/download_%J.out
#SBATCH --partition=amilan

module load sra-toolkit

mkdir -p logs

SRR_val=SRR10088651
sample_name=control_5

scratch_path=/scratch/alpine/$USER/analysis/wells/analysis/sussel/rbfox2/rbfox2_analysis/rna_seq_reanalysis/raw_data
mkdir -p $scratch_path

wget https://sra-pub-run-odp.s3.amazonaws.com/sra/${SRR_val}/${SRR_val}
fastq-dump -I --split-files --gzip --origfmt ${SRR_val}
mv ${SRR_val}_1.fastq.gz 	$scratch_path/${sample_name}_R1.fastq.gz
mv ${SRR_val}_2.fastq.gz 	$scratch_path/${sample_name}_R2.fastq.gz
rm ${SRR_val}

ln -s $scratch_path/${sample_name}_R1.fastq.gz ./
ln -s $scratch_path/${sample_name}_R2.fastq.gz ./