library(tidyverse)
library(here)

args <- commandArgs(trailingOnly = TRUE)

results_dir <- args[1]
html_file <- args[2]

quality_rmd <- here("src/scripts/quality_analysis_plots.Rmd")

# Save files
mapping_res <- here("results/mapping_quality.csv")
dup_file <- here("results/all_quality_res.csv")

# Read in results
bowtie_dir <- file.path(results_dir, "bowtie2_cutadapt_trim")

sample_info <- read.table(here("files/samples.tsv"), sep = "\t", header = TRUE)
samples <- sample_info$sample

bowtie_res <- read.csv(file.path(bowtie_dir, "bowtie2_stats.csv"))
dup_res <- read.table(file.path(bowtie_dir, "dup_stats.txt"), sep = "\t",
                      header = TRUE)

mito_res <- lapply(samples, function(sample){
  mito_single <- read.csv(file.path(bowtie_dir,
                                    paste0(sample, "_mito_stats.txt")))
  return(mito_single)
})

mito_res <- do.call(rbind, mito_res)

read_counts <- bowtie_res %>%
  dplyr::select(total_reads, sample) %>%
  mutate(label_reads = round(total_reads / 1000000)) %>%
  mutate(label_reads = str_c(label_reads, " million")) %>%
  dplyr::distinct()

# Fix alignments to give more information
bowtie_res$concordant_alignment <- (bowtie_res$one_alignment + 
                                      bowtie_res$multi_alignment)

bowtie_res$overall_alignment <- as.numeric(gsub("%", "",
                                                bowtie_res$overall_alignment))

bowtie_res$total_reads_aligned <- (bowtie_res$overall_alignment / 100 *
                                     bowtie_res$total_reads)

bowtie_res$non_concordant_alignment <- as.character((bowtie_res$total_reads_aligned - 
                                              bowtie_res$concordant_alignment) /
                                             bowtie_res$total_reads * 100)

bowtie_res$no_alignment_percent <- as.character(100 - bowtie_res$overall_alignment)
  
all_results <- bowtie_res %>%
  pivot_longer(cols = c(no_alignment_percent, 
                        one_alignment_percent, multi_alignment_percent,
                        non_concordant_alignment),
               names_to = "alignment_type", values_to = "read_percent") %>%
  dplyr::mutate(read_percent_fixed = gsub("\\(", "", read_percent)) %>%
  dplyr::mutate(read_percent_fixed = gsub("%\\)", "", read_percent_fixed)) %>%
  dplyr::mutate(read_percent = as.numeric(read_percent_fixed))


write.csv(all_results, mapping_res)

# Clean up mito and duplication removal
full_results <- merge(bowtie_res, dup_res, by = "sample")
full_results <- merge(mito_res, full_results, by = "sample")

# In full results
# Mito mapped reads / 2 = concordant alignment 
# Ignore all reads
full_results$percent_mito <- 100 - (full_results$no_mito_reads /
                                      full_results$mapped_reads * 100)


write.csv(full_results, dup_file)

# Read in fastq screen as well
custom_params <- list(
  sample_table = here("files/samples.tsv"),
  mapping_res = mapping_res,
  fastqc1 = file.path(results_dir, "fastqc_pre_trim_summary_untrimmed.tsv"),
  fastqc2 = file.path(results_dir, "fastqc_cutadapt_summary_trimmed.tsv"),
  duplicates = dup_file,
  fastq_screen = file.path(results_dir, "fastq_screen", 
                           "fastq_screen_info.csv"),
  snakemake = TRUE
)

rmarkdown::render(quality_rmd, output_file = html_file, params = custom_params,
                  envir = new.env())
                
