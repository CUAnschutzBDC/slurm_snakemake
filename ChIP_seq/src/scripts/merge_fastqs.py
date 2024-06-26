from collections import defaultdict
import argparse
import pandas as pd
import os
import subprocess
import gzip

def main():
	opts = setup()
	r1_dict, r2_dict = build_dictionary(opts.sample_sheet)
	merge_files(r1_dict, opts.out_dir, "R1")
	merge_files(r2_dict, opts.out_dir, "R2")

#########
# Setup #
#########

def setup():
	"""
	Gets command line arguments and returns a Namespace object
	"""

	# House keeping to read in arguments from the command line
	parser = argparse.ArgumentParser()

	parser.add_argument("-s", "--sample_sheet", dest = "sample_sheet",
		help = "The sample sheet passed to the snakefile",
		action = "store",
		metavar = "\b")

	parser.add_argument("-o", "--out", dest = "out_dir",
		help = "The directory to put the output",
		action = "store",
		metavar = "\b")

	args = parser.parse_args()

	return(args)

def build_dictionary(sample_table):
	r1_dict = defaultdict(list)
	r2_dict = defaultdict(list)
	sample_info = pd.read_table(sample_table).set_index("sample", drop = False)
	for index, row in sample_info.iterrows():
		fastq1 = os.path.join(row["data_dir"], row["fastq1"])
		fastq2 = os.path.join(row["data_dir"], row["fastq2"])
		r1_dict[row["sample"]].append(fastq1)
		r2_dict[row["sample"]].append(fastq2)

	return(r1_dict, r2_dict)

def merge_files(read_dict, output_dir, read):
	os.makedirs(output_dir, exist_ok = True)

	for sample in read_dict:
		fastq_file =  os.path.join(output_dir, f"{sample}_{read}.fastq.gz")
		if len(read_dict[sample]) > 1:
			with gzip.open(fastq_file, "wb") as fastq:
				for file_name in read_dict[sample]:
					with gzip.open(file_name, "rb") as in_fastq:
						for line in in_fastq:
							fastq.write(line)
		else:
			subprocess.run(["ln", "-s", read_dict[sample][0], fastq_file], stderr=subprocess.PIPE, check=True)			


if __name__ == "__main__":
	main()