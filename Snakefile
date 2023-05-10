""" Snake pipeline for running RNAseq analysis """
results = "results"
samples = ["s1", "s2"]
# Final output files
rule all:
    input:
        # Determine what needs to be run by copying completed files
        expand(
            "{results}/{sample}_two.txt",
            results = results, sample = samples
            )

        
# Snakes to run
include: "src/rules/test_rules.snake"
