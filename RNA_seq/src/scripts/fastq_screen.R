library(tidyverse)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)

save_dir <- args[1]
file_pattern <- "_S[0-9]+_L[0-9]*_R[1|2]_[0-9]*_screen.txt"
#file_pattern <- "_R[1|2]_screen.txt"
all_files <-list.files(path = save_dir, pattern = file_pattern, full.names = FALSE)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

all_results <- lapply(all_files, function(file_name){
	#file_name <- "Input_1_S15_L006_R1_001_screen.txt"
	sample <- gsub(file_pattern, "", file_name)

	file <- read.table(file.path(save_dir, file_name), sep = "\t",
					fill = TRUE, comment.char = "", header = TRUE, 
					skip = 1)

	# Remove the final line which is unmapped
	main_file <- file[1:(nrow(file) -1),]
	genome <- main_file[ , "Genome"]
	main_file$Genome <- NULL
	no_hits <- file[(nrow(file)):nrow(file), 1]
	no_hits <- strsplit(no_hits, ": ")[[1]][2]

	# Make numeric
	main_file[] <- lapply(main_file, as.numeric)
	main_file$Genome <- genome

	# Want 1 barplot for unmapped
	# Want 1 barplot for hits
	# Here we have all 1 hit 1 genome, all multiple hits one genome, combined one hit multiple
	# genomes and multiple hits multiple genomes
	main_file$percent <- main_file$X.One_hit_one_genome.1 + main_file$X.Multiple_hits_one_genome.1
	main_file <- main_file[ , c("Genome", "percent")]
	main_file <- rbind(main_file, data.frame("Genome" = "no_hits", "percent" = as.numeric(no_hits)))

	# The multi mappers are hard to determine because there are so many overlaps in the column. 
	# Here, I just assume that any reads not accounted for are multimapped
	main_file <- rbind(main_file, data.frame("Genome" = "multiple_genomes", "percent" =100 - sum(main_file$percent)))

	main_file$sample <- sample

	return(main_file)

})

# Combine all files
all_results <- do.call(rbind, all_results)

# Make the plotting levels better?
all_results$Genome <- factor(all_results$Genome, levels = unique(all_results$Genome))

# Pick some colors
colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(name = "Set1", n = 9))(length(levels(all_results$Genome)))

names(colors) <- levels(all_results$Genome)

# Make a barplot
barplot <- ggplot2::ggplot(all_results, ggplot2::aes(x = sample, y = percent, fill = Genome)) +
	ggplot2::geom_bar(stat = "identity", position = "stack") +
	ggplot2::coord_flip() +
	ggplot2::scale_fill_manual(values = colors)

fig_size <- length(all_files) / 1.25

write.csv(all_results, file.path(save_dir, "fastq_screen_info.csv"))

# Make the file
pdf(file.path(save_dir, "fastq_screen_barplot.pdf"),
    width = 8, height = fig_size)
print(barplot)
dev.off()
