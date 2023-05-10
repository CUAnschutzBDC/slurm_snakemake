# Slurm snakemake

This repository contains the code necessary to run snakemake in a slurm environment. It includes the profiles necessary to run snakemake with or without singularity.

## Usage

To use this repo, first clone it. Once you have cloned it, update the Snakefile and rules to fit your pipeline. To then run the snakemake pipeline either run the snakecharmer:

```bash
sbatch snakecharmer.sh
```

or you can run directly from the command line using

```bash
snakemake --profile {path/to/profile)
```

Replace the `{path/to/profile}` with either `slurm-profile/` or `slurm-profile-singularity` depending on if you are using singularity.

## Planned additions
Currently, submitting with the profile fails when trying to run groups because snakemake cannot identify the needed resources. More testing and reading is required to figure this out.
