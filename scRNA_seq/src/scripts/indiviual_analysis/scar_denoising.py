# Follows https://github.com/KenLauLab/dropkick/blob/master/dropkick_tutorial.ipynb

import scanpy as sc; sc.set_figure_params(color_map="viridis", frameon=False)
import scar
import os
import matplotlib.pyplot as plt

def main():
	results_dir = "/beevol/home/wellskri/Analysis/Mia_Smith/Catherine_Nicolas/full_antigen_pos_data/results"
	sample = "JH310-12_AP"

	results_dir = snakemake.wildcards["results"]
	sample = snakemake.wildcards["sample"]
	raw_data_path = os.path.join(results_dir, sample, "outs", "multi", "count", "raw_feature_bc_matrix.h5") 
	filtered_data_path = os.path.join(results_dir, sample, "outs", "per_sample_outs", sample, "count",
		                              "sample_filtered_feature_bc_matrix.h5")
	save_dir = os.path.join(results_dir, "R_analysis", sample)
	image_dir = os.path.join(results_dir, "R_analysis", sample, "images", "scar")
	file_dir = os.path.join(results_dir, "R_analysis", sample, "files")

	adata = sc.read_10x_mtx(data_path, var_names = "gene_symbols")

	# simple preprocessing of anndata object to get metrics for plot
	adata = dk.recipe_dropkick(adata, n_hvgs=None, X_final="raw_counts")

	# Make QC plots
	qc_plt = dk.qc_summary(adata)

	plt.savefig(os.path.join(image_dir, "qc_plot.pdf"), format='pdf')

	# Run dropkick pipeline function
	adata_model = dk.dropkick(adata, n_jobs=5)


	score_plt = dk.score_plot(adata)

	plt.savefig(os.path.join(image_dir, "score_plot.pdf"), format='pdf')


	coef_plt = dk.coef_plot(adata)
	plt.savefig(os.path.join(image_dir, "coef_plot.pdf"), format='pdf')

	# Filter to only cells considered "real"
	adata = adata[(adata.obs.dropkick_label == "True"),:]

	# Save the file
	adata.obs.to_csv(os.path.join(file_dir, "dropkick_cells.csv"))

if __name__ == "__main__":
	main()