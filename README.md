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

### Groups
Running snakemake with groups will allow us to submit fewer jobs to the cluster. For instance, in a bulk RNA-seq pipeline, adaptor trimming, mapping, and feature counting for one sample does not need to be separate jobs because they can't be run in parallel. Here, they can all be added to a group and submitted only one time per sample. An example of running this way is in `src/rules/test_rules_resources.snake`. This is included below

```bash
rule one_rule:
	output:
		"{results}/{sample}_one.txt"
	params:
		tmp_dir = "/tmp"
	resources:
		job_name="test",
		mem=1000,
		time="0:10:00",
		partition="amilan",
		slurm=lambda wildcards: f"output={wildcards.results}/logs/group1/{wildcards.sample}_test.out error={wildcards.results}/logs/group1/{wildcards.sample}_test.err"
	threads:
		1
	group:
		"group1"
	shell:
		"""
		echo "starting one_rule"
		sleep 15
		touch {output}
		echo "finished one_rule"
		"""


rule two_rule:
	input:
		"{results}/{sample}_one.txt"
	output:
		"{results}/{sample}_two.txt"
	params:
		tmp_dir = "/tmp"
	resources:
		job_name="test",
		mem=1000,
		time="0:10:00",
		partition="amilan",
		slurm=lambda wildcards: f"output={wildcards.results}/logs/group1/{wildcards.sample}_test.out error={wildcards.results}/logs/group1/{wildcards.sample}_test.err"
	threads:
		1
	group:
		"group1"
	shell:
		"""
		echo "starting two_rule"
		sleep 15
		touch {output}
		echo "finished two_rule"
		"""
```

A few notes on using resources:

* All rules grouped together must have the same "group" parameter. Resources in a group must contain the same name and partition. If the time is a string as in this example, than that must be the same for all members of the group as well.
* To add log files, you can pass these to slurm with `slurm="slurm-arg={arg}"` For example the log file is `slurm="output=path/to/log.out". Some slurm arguments are automatically included while some are not. More description under ["Rule specific resource configuration"](https://github.com/Snakemake-Profiles/slurm). This log file must be identical for all grouped rules.
* To add any wildcard information to the `resources` `lambda` must be used: `job_name=lambda wildcards: f"{wildcards.sample}-rule"`. An example of this is seen for the log files in the example above.
