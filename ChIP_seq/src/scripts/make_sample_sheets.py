import os
import argparse
import re

def main():
	opts = setup()

	# make the sample_info file
	sample_list = make_sample_sheet(opts.in_dir, opts.out_dir)
	
#########
# Setup #
#########

def setup():
    """
    Gets command line arguments and returns a Namespace object
    """

    # House keeping to read in arguments from the command line
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", dest = "in_dir",
        help = "The directory that holds the log files, default is raw_data",
        default = "raw_data",
        action = "store",
        metavar = "\b")

    parser.add_argument("-o", "--out", dest = "out_dir",
        help = "The directory to put the output, default is files",
        default = "files",
        action = "store",
        metavar = "\b")

    args = parser.parse_args()

    return(args)

##############
# Make files #
##############
def make_sample_sheet(in_dir, out_dir):
	"""
	Takes all files in an input directory and finds R1 files by searching for R1_001 in the file name
	R2 file and sample are determined based on the R1 file. This assumes the normal fastq structure output
	by illumina bcl2fastq

	Also builds a list of samples so the sample_info and rmats_sample_info can be built
	"""
	# Generate an empty sample list
	sample_list = []
	out_file = os.path.join(out_dir, "samples.tsv")

	# First list all fastq files in the directory
	files = os.listdir(in_dir)
	# Keep only R1 files
	filtered_files = [file_name for file_name in files if not file_name.startswith("._") and "_R1_001" in file_name]

	# Set fastq pattern
	fastq_pattern = r"_S[0-9]*_L[0-9]*_R[1|2]_001\.fastq\.gz"
	with open(out_file, "w") as output:
		output.write("sample\tfastq1\tfastq2\tdata_dir\tcontrol\treplicate\n")
		for R1_file in filtered_files:

			# Determine R2 file and sample based on the R1 file
			R2_file = re.sub("_R1_001", "_R2_001", R1_file)
			sample = re.sub(fastq_pattern, "", R1_file)

			# Add sample to sample list
			sample_list.append(sample)

			# Make list of items to write to file
			write_list = [sample, R1_file, R2_file, in_dir]

			# Write all info to the output file
			output.write("\t".join(write_list) + "\t\t\n")
			print("Warning: samples requires additional fields to be filled before running!")
			print("Include what control to pair for each sample or 'No_ctl' if it is the control")
			print("Include the replicate information")

	return sample_list


if __name__ == "__main__":
	main()