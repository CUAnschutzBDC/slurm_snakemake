library(tidyverse)
library(here)

args <- commandArgs(trailingOnly = TRUE)

results_dir <- args[1]
project <- args[2]
html_file <- args[3]

quality_rmd <- here("src/scripts/QC_only.Rmd")

custom_params <- list(
  genome = "GRCm38",
  DE_alpha = 0.05,
  DE_lfc = 1,
  counts = here(results_dir, paste0(project, "_countTable_cutadapt_trim.txt")),
  fastqc_one = here(results_dir, "fastqc_pre_trim_summary_untrimmed.tsv"),
  fastqc_two = here(results_dir, "fastqc_cutadapt_summary_trimmed.tsv"),
  star_stats = here(results_dir, "star_summary_cutadapt_trim.tsv"),
  output_dir = here(results_dir, "R_analysis_check"),
  sample_info = here("files", "sample_info.csv"),
  fastq_screen = here(results_dir, "fastq_screen", "fastq_screen_info.csv"),
  sample_column = "genotype",
  snakemake = TRUE
)

rmarkdown::render(quality_rmd, output_file = html_file, params = custom_params,
                  envir = new.env())
