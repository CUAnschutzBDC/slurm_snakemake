library(here)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

in_dir <- args[1]
output_1 <- args[2]
output_2 <- args[3]

all_files <- list.files(in_dir)

merge_files <- function(all_files, in_dir, type, output_file){
	merge_files <- all_files[grepl(paste0("MATS.", type, ".txt"), all_files)]
	all_files <- lapply(merge_files, function(file_name){
		splice_type <- gsub("\\..*", "", file_name)
		read_file <- read.table(file.path(in_dir, file_name), header = 1)
		read_file$splice_type <- splice_type
		return(read_file)
	})

	all_files <- dplyr::bind_rows(all_files)
	write.table(all_files, file = output_file)
}

merge_files(all_files, type = "JCEC", in_dir = in_dir, output_file = output_1)
merge_files(all_files, type = "JC", in_dir = in_dir, output_file = output_2)