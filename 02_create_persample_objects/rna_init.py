import os
import argparse
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
import seaborn as sns
import anndata
from scipy.io import mmread
import mymodule.utils as utils
from matplotlib.backends.backend_pdf import PdfPages

def check_files_exist(sample_name, outs_dir):
    required_files = [
        os.path.join(outs_dir, sample_name, "plots", "pre_filt_QCplots.pdf"),
        os.path.join(outs_dir, sample_name, "plots", "cr_filt_QCplots.pdf"),
        os.path.join(outs_dir, sample_name, "pre_filt_rna_adata.h5ad"),
        os.path.join(outs_dir, sample_name, "cr_filt_rna_adata.h5ad")
    ]
    return all([os.path.exists(file_path) for file_path in required_files])

def ribo_content_plots(adata):

    fig, axes = plt.subplots(1, 3, figsize=(8*3, 6))

    axes[0].hist(adata.obs['ribo_pct'], bins=100)
    axes[0].set_xlabel('% Ribo-content', fontsize=14)
    axes[0].set_ylabel('Frequency', fontsize=14)

    axes[1].scatter(adata.obs['log1p_total_counts'], adata.obs['ribo_pct'])
    axes[1].set_xlabel('Log library size', fontsize=14)
    axes[1].set_ylabel('% Ribo-content', fontsize=14)

    axes[2].scatter(adata.obs['log1p_n_genes_by_counts'], adata.obs['ribo_pct'])
    axes[2].set_xlabel('Log num. genes per cell', fontsize=14)
    axes[2].set_ylabel('% Ribo-content', fontsize=14)

    plt.tight_layout()
    return fig

def mt_content_plots(adata):

    fig, axes = plt.subplots(1, 3, figsize=(8*3, 6))

    axes[0].hist(adata.obs['mito_pct'], bins=100)
    axes[0].set_xlabel('% MT-content', fontsize=14)
    axes[0].set_ylabel('Frequency', fontsize=14)

    axes[1].scatter(adata.obs['log1p_total_counts'], adata.obs['mito_pct'])
    axes[1].set_xlabel('Log library size', fontsize=14)
    axes[1].set_ylabel('% MT-content', fontsize=14)

    axes[2].scatter(adata.obs['log1p_n_genes_by_counts'], adata.obs['mito_pct'])
    axes[2].set_xlabel('Log num. genes per cell', fontsize=14)
    axes[2].set_ylabel('% MT-content', fontsize=14)

    plt.tight_layout()
    return fig

def genes_metadata_plots(n_cells_by_counts):
    fig, axes = plt.subplots(1, 3, figsize=(8*3, 6))
    
    axes[0].hist(n_cells_by_counts, bins=100)
    axes[0].set_xlabel('Number of cells a gene is expressed in', fontsize=14)
    axes[0].set_ylabel('Frequency', fontsize=14)
    axes[0].set_title('Histogram of number of cells each gene is expressed in', fontsize=14)

    axes[1].hist(np.log(n_cells_by_counts + 1), bins=100)
    axes[1].set_xlabel('Log - Number of cells a gene is expressed in', fontsize=14)
    axes[1].set_ylabel('Frequency', fontsize=14)
    axes[1].set_title('Histogram of log of number of cells each gene is expressed in', fontsize=14)

    axes[2].hist(np.log(n_cells_by_counts + 1), bins=100)
    axes[2].set_xlabel('Log - Number of cells a gene is expressed in', fontsize=14)
    axes[2].set_ylabel('Frequency, ylim trimmed', fontsize=14)
    axes[2].set_title('Histogram of log of number of cells each gene is expressed in', fontsize=14)
    axes[2].set_ylim([0, 1000])

    plt.tight_layout()
    return fig

def cells_genes_correlation_plots(log1p_total_counts, log1p_n_genes_by_counts):
    fig, axes = plt.subplots(1, 3, figsize=(8*3, 6))
    
    sns.histplot(log1p_total_counts, bins=100, ax=axes[0])
    axes[0].set_xlabel('Log-library size', fontsize=14)
    axes[0].set_ylabel('Frequency', fontsize=14)
    axes[0].set_title('Histogram of log library size', fontsize=14)

    sns.histplot(log1p_n_genes_by_counts, bins=100, ax=axes[1])
    axes[1].set_xlabel('Log of num. genes per cell', fontsize=14)
    axes[1].set_ylabel('Frequency', fontsize=14)
    axes[1].set_title('Histogram of number of genes expressed in each cell', fontsize=14)

    sns.scatterplot(x=log1p_total_counts, y=log1p_n_genes_by_counts, ax=axes[2])
    axes[2].set_ylabel('Log of num. genes per cell', fontsize=14)
    axes[2].set_xlabel('Log library size', fontsize=14)
    corr_coef = np.corrcoef(log1p_total_counts, log1p_n_genes_by_counts)[0, 1]
    axes[2].set_title(f'Correlation = {corr_coef:.3f}', fontsize=14)

    plt.tight_layout()
    return fig

def load_anndata(mtx_dir):
    # load in mtx, barcodes and features
    matrix_path = os.path.join(mtx_dir, 'matrix.mtx.gz')
    matrix = mmread(matrix_path).tocsr()
    barcodes_path = os.path.join(mtx_dir, 'barcodes.tsv.gz')
    barcodes = pd.read_csv(barcodes_path,sep='\t',compression='gzip',header=None,
                          names = ['og_barcode'])
    features_path = os.path.join(mtx_dir, 'features.tsv.gz')
    features = pd.read_csv(features_path, sep='\t',compression='gzip',header=None, 
                           names=['ensemblID','name', 'feature_type', 'chr', 'start', 'end'])
    
    # select gene features
    gene_features = np.where(features['feature_type'] == 'Gene Expression')[0]
    
    # generate anndata
    var = features.iloc[gene_features]
    var.index = var['name'].tolist()
    obs = barcodes
    obs.index = obs['og_barcode'].tolist()
    adata = anndata.AnnData(matrix[gene_features,:].T, var = var, obs = obs)
    adata.var_names_make_unique()
    return adata

def preprocess_anndata(sample_name, cr_outs, outs_dir, overwrite=False):
    
    print('Processing sample', sample_name)
    mtx_dir = os.path.join(cr_outs, "raw_feature_bc_matrix")

    # Check if all required files are already present
    if os.path.exists(mtx_dir):
        adata = load_anndata(mtx_dir)
        # rename barcodes to sample_name+barcodes and add the sample data to anndata
        adata.obs_names = [f"{sample_name}#{barcode}" for barcode in adata.obs_names]
        adata.obs['batch'] = sample_name

        # get qc metrics
        adata = utils.calculate_qc_metrics(adata)
        # save as a new h5ad file
        sample_outs_dir = os.path.join(outs_dir, sample_name)
        os.makedirs(sample_outs_dir, exist_ok=True)
        adata.write_h5ad(os.path.join(sample_outs_dir, f"pre_filt_rna_adata.h5ad"))
        print("Saving modified adata as pre_filt_rna_adata.h5ad", sample_name)

        # create a new directory to save the plots
        plots_dir = os.path.join(sample_outs_dir, "plots")
        os.makedirs(plots_dir, exist_ok=True)

        # generate qc plots prior to filteration
        print('Generating pre-filtering QC plots', sample_name)
        pdf_file = os.path.join(plots_dir, 'pre_filt_QCplots.pdf')
        with PdfPages(pdf_file) as pdf:
            pdf.savefig(cells_genes_correlation_plots(adata.obs['log1p_total_counts'], adata.obs['log1p_n_genes_by_counts']))
            pdf.savefig(genes_metadata_plots(adata.var['n_cells_by_counts']))
            pdf.savefig(mt_content_plots(adata))
            pdf.savefig(ribo_content_plots(adata))
                
        # get cellranger filtered barcodes
        print('Subsetting by cellranger filtered cells', sample_name)
        filt_bc_path = os.path.join(cr_outs, "filtered_feature_bc_matrix/barcodes.tsv.gz")
        cr_cells = pd.read_csv(filt_bc_path,sep='\t',compression='gzip',header=None)[0].tolist()
        cr_cells = [sample_name + '#' + x for x in cr_cells]
                
        # subset + embed
        filt_ad = utils.subset_and_reprocess_rna(adata,cr_cells,remove_mito_ribo = False,
                             filter_genes = False)
                
        # run doubletdetection on filtered anndata
        utils.run_doubletdetection(filt_ad, sample_col = 'batch', layer = 'raw_counts')
        
        # generate qc plots prior to filteration
        print('Generating cellranger-filtered QC plots', sample_name)
        pdf_file = os.path.join(plots_dir, 'cr_filt_QCplots.pdf')
        with PdfPages(pdf_file) as pdf:
            pdf.savefig(cells_genes_correlation_plots(filt_ad.obs['log1p_total_counts'], filt_ad.obs['log1p_n_genes_by_counts']))
            pdf.savefig(genes_metadata_plots(filt_ad.var['n_cells_by_counts']))
            pdf.savefig(mt_content_plots(filt_ad))
            pdf.savefig(ribo_content_plots(filt_ad))
                
        # save filtered anndata
        filt_ad.write_h5ad(os.path.join(sample_outs_dir, f"cr_filt_rna_adata.h5ad"))
        print('Finished processing sample', sample_name)
        print('Output files are stored in:', sample_outs_dir)

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process data for quality control.')
    parser.add_argument('--sample', '-s', type=str, help='Sample name (should match with the directory name)')
    parser.add_argument('--data', '-i', type=str, help='Path to cellranger output directory.')
    parser.add_argument('--out', '-o', type=str, help='Path to the output directory.')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite existing files.')
    args = parser.parse_args()
    sample_name = args.sample
    cr_outs = args.data
    outs_dir = args.out
    overwrite = args.overwrite
    
    if overwrite == False:
        exists = check_files_exist(sample_name, outs_dir)
        if exists == True:
            print("Results for sample", sample_name, "is already present. If you would like to overwrite the results please add overwrite option while running the script.")
        else:
            preprocess_anndata(sample_name, cr_outs, outs_dir)
    else:
        preprocess_anndata(sample_name, cr_outs, outs_dir)