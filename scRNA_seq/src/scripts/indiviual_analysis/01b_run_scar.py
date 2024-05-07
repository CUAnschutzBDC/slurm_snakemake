import os
import scanpy as sc
import scar
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import argparse

def main():
	options = setup()

	results_dir = options.results_dir
	sample = options.sample_name
	raw_data_path = os.path.join(results_dir, sample, "outs", "per_sample_outs", 
		                         sample, "count", "raw_feature_bc_matrix.h5") 

	filtered_data_path = os.path.join(results_dir, sample, "outs", "per_sample_outs",
	                                  sample, "count",
	                                  "sample_filtered_feature_bc_matrix.h5")
	save_dir = os.path.join(results_dir, "R_analysis", sample)
	image_dir = os.path.join(results_dir, "R_analysis", sample, "images", "scar")
	file_dir = os.path.join(results_dir, "R_analysis", sample, "files")

	save_file = os.path.join(results_dir, "R_analysis", sample, "files",
	                         "scar_denoised.csv")

	sample_raw = sc.read_10x_h5(raw_data_path, gex_only = False)
	sample_raw.var_names_make_unique()

	sample_filtered = sc.read_10x_h5(filtered_data_path, gex_only = False)
	sample_filtered.var_names_make_unique()

	scar.setup_anndata(
	    adata = sample_filtered,
	    raw_adata = sample_raw,
	    prob = 0.7, # try a lower prob value when error ocurs
	    min_raw_counts=2,
	    kneeplot = True
	)

	ambient_profile = sample_filtered.uns["ambient_profile_Antibody Capture"]




	warnings.simplefilter("ignore")

	# Need to pull out raw ADT counts
	raw_obj = sample_filtered.X
	dense_array = raw_obj.toarray()

	# Convert the dense array to a DataFrame
	dense_df = pd.DataFrame(dense_array)
	dense_df.index = sample_filtered.obs.index
	dense_df.columns = sample_filtered.var.index

	# Pull out  ADTs only 
	gene_info = sample_filtered.var
	gene_info = gene_info[gene_info["feature_types"] == "Antibody Capture"]

	subset_df = dense_df[gene_info.index]


	ADT_scar = scar.model(raw_count = subset_df,
	                 ambient_profile = ambient_profile,  # Providing ambient profile is recommended for CITEseq; in other modes, you can leave this argument as the default value -- None
	                 feature_type = 'ADT', # "ADT" or "ADTs" for denoising protein counts in CITE-seq
	                 count_model = 'binomial',   # Depending on your data's sparsity, you can choose between 'binomial', 'possion', and 'zeroinflatedpossion'
	                )

	ADT_scar.train(epochs=500,
	               batch_size=64,
	               verbose=True
	              )

	# After training, we can infer the true protein signal
	ADT_scar.inference()  # by defaut, batch_size=None, set a batch_size if getting GPU memory issue

	denoised_ADT = pd.DataFrame(ADT_scar.native_counts, index=subset_df.index, columns=subset_df.columns)

	denoised_ADT.to_csv(save_file)


#########
# Setup #
#########

def setup():
	"""
	Gets command line arguments and returns a Namespace object
	"""

	# House keeping to read in arguments from the command line
	parser = argparse.ArgumentParser()

	parser.add_argument("-d", "--directory", dest = "results_dir",
		help = "Path to the output directory used by snakemake",
		default = "none",
		action = "store",
		metavar = "\b")

	parser.add_argument("-s", "--sample_name", dest = "sample_name",
		help = "The name of the sample",
		default = "none",
		action = "store",
		metavar = "\b")

	args = parser.parse_args()

	return(args)


if __name__ == "__main__":
	main()
