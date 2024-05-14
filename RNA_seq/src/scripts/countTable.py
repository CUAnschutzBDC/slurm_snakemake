from collections import defaultdict
import re
import os

file_list = snakemake.params["in_files"]
output_file = snakemake.output[0]

trim_method = snakemake.wildcards.trim_method
results = snakemake.wildcards.results

gene_dict = defaultdict(list)
sample_list = list()

#file_path = results + "/featureCount" + "_" + trim_method + "_trim/"

for i in file_list:
	sample_name = os.path.split(i)[-1]
	sample_name = re.sub(r'_countsOutput', '', sample_name)
	sample_list.append(sample_name)
	with open(i, "r") as countFile:
		for line in countFile:
			line = line.strip().split("\t")
			if "ENS" in line[0]:
				ens_id = line[0]
				gene = line[6]
				count = line[8]
				gene_ens = gene + "_" + ens_id
				gene_dict[gene_ens].append(count)

with open(output_file, "w") as count_file:
	count_file.write("gene" + "\t" + "\t".join(sample_list) + "\n")
	for gene in gene_dict:
		count_file.write(gene + "\t" + "\t".join(gene_dict[gene]) + "\n")
