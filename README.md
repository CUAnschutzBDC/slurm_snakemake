# Analysis pipelines
Writen by Kristen Wells and Chris Hill

A collection of snakemake pipelines to analyze omics datasets. Current pipelines include:

* A scRNA-sequencing pipeline that analyzed 10x genomics datasets. Can be used with scRNA-seq, scCITE-seq, scVDJ-seq, and hashtagging (and any combination of those).
* A basic bulk RNA-sequencing pipeline that includes options for trimming (cutadapt, bbduk or none). The pipeline includes alignment with star and read counting with FeatureCounts. This pipeline will output a quality Rmarkdown as the final step.
* A pipeline to process ChIP-seq data that includes options for trimming (cutadapt, bbduk, or none). The pipeline includes alignment iwth bowtie2, peak calling with Macs2 and initial analysis with diffbind. This pipeline will output a quality Rmarkdown as the final step.

Pipelines that will be added soon include:
* An ATAC-seq pipeline

Snakemake pipelines here were written to work on a slurm cluster but should be easily ported to other systems. All pipelines rely on docker and singularity images.

Snakemake must be installed to use these pipelines.

1. Download and install miniconda3: For Linux
```{bash}
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh bash Miniconda3-latest-Linux-x86_64.sh
```
2. Install Snakemake:
```{bash}
conda install snakemake -c bioconda -c conda-forge
```

Steps to run indivdual pipelines are within their own directories.

If you use any of the pipelines in this repository, please give appropriate credit.
