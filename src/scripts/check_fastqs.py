#!/usr/bin/env python3
from Bio import SeqIO
import sys
import argparse
from pathlib import Path
import gzip

def main():
	args = setup()

	try:
		reads_match, read_count = validate_paired_fastq(args.r1, args.r2, args.verbose)
		
		print(f"\nProcessed {read_count:,} reads")
		if reads_match:
			print("All read names match between R1 and R2 files")
			sys.exit(0)
		else:
			print("Found mismatches between R1 and R2 read names")
			sys.exit(1)
			
	except Exception as e:
		print(f"Error: {str(e)}")
		sys.exit(1)


def validate_paired_fastq(r1_path: str, r2_path: str, verbose: bool = False) -> tuple[bool, int]:
	"""
	Validate that read names match between paired FASTQ files.

	Args:
		r1_path: Path to R1 FASTQ file
		r2_path: Path to R2 FASTQ file
		verbose: If True, print mismatched read names
		
	Returns:
		tuple containing:
			- Boolean indicating if all reads match
			- Number of reads processed
	"""
	# Check if files exist
	if not Path(r1_path).exists():
		raise FileNotFoundError(f"R1 file not found: {r1_path}")
	if not Path(r2_path).exists():
		raise FileNotFoundError(f"R2 file not found: {r2_path}")

	# Open both files using SeqIO
	with gzip.open(r1_path, "rt") as f1, gzip.open(r2_path, "rt") as f2:
		r1_parser = SeqIO.parse(f1, "fastq")
		r2_parser = SeqIO.parse(f2, "fastq")
		
		reads_match = True
		read_count = 0
		
		# Iterate through both files simultaneously
		for r1_read, r2_read in zip(r1_parser, r2_parser):
			read_count += 1
			if read_count % 1000000 == 0 and verbose:
				print(str(read_count) + " reads processed")
			
			# Get read names (IDs)
			r1_name = r1_read.id
			r2_name = r2_read.id
			
			# Check if names match
			if r1_name != r2_name:
				reads_match = False
				if verbose:
					print(f"Mismatch at read {read_count}:")
					print(f"R1: {r1_name}")
					print(f"R2: {r2_name}")
					print("---")
					
				sys.exit()
				
		# Check if one file has more reads than the other
		try:
			next(r1_parser)
			print("Warning: R1 file contains more reads than R2")
			reads_match = False
		except StopIteration:
			try:
				next(r2_parser)
				print("Warning: R2 file contains more reads than R1")
				reads_match = False
			except StopIteration:
				pass

	return reads_match, read_count

#########
# Setup #
#########

def setup():
	"""
	Gets command line arguments and returns a Namespace object
	"""

	# House keeping to read in arguments from the command line
	parser = argparse.ArgumentParser(description="Validate paired FASTQ files")
	parser.add_argument("-r1", help="Path to R1 FASTQ file")
	parser.add_argument("-r2", help="Path to R2 FASTQ file")
	parser.add_argument("-v", "--verbose", action="store_true", 
						help="Print mismatched read names")
	args = parser.parse_args()

	return(args)

if __name__ == "__main__":
    main()