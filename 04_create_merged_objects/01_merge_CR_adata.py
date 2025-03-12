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
from pathlib import Path


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

def load_anndatas(sample_names, input_dir):
    # load both unfiltered and filtered anndatas
    anndata_list = []
    for sample in sample_names:
        mtx_dir = os.path.join(input_dir,sample,"filtered_feature_bc_matrix")
        adata = load_anndata(mtx_dir)
        adata.obs['og_ident'] = sample
        adata.obs.index = adata.obs['og_ident']+'#'+adata.obs['og_barcode']
        #rna_ad_path = os.path.join(input_dir,sample,'cr_filt_rna_adata.h5ad')
        anndata_list.append(adata)
    return anndata_list

def combine_RNA_anndata(anndata_list):
    
    # prep raw counts, obs and var
    count_mtx = scipy.sparse.vstack([ad.X for ad in anndata_list])
    obs_df = pd.concat([ad.obs for ad in anndata_list])
    var_df = anndata_list[0].var
       
    # join together
    comb_ad = anndata.AnnData(count_mtx, 
                          obs = obs_df, 
                          var = var_df)
    sc.pp.calculate_qc_metrics(comb_ad,inplace=True)
    return comb_ad

def kneepoint(vec):
    #computes kneepoint in PCA cumulative variance explained
    curve =  [1-x for x in vec]
    nPoints = len(curve)
    allCoord = np.vstack((range(nPoints), curve)).T
    np.array([range(nPoints), curve])
    firstPoint = allCoord[0]
    lineVec = allCoord[-1] - allCoord[0]
    lineVecNorm = lineVec / np.sqrt(np.sum(lineVec**2))
    vecFromFirst = allCoord - firstPoint
    scalarProduct = np.sum(vecFromFirst * np.tile(lineVecNorm, (nPoints, 1)), axis=1)
    vecFromFirstParallel = np.outer(scalarProduct, lineVecNorm)
    vecToLine = vecFromFirst - vecFromFirstParallel
    distToLine = np.sqrt(np.sum(vecToLine ** 2, axis=1))
    idxOfBestPoint = np.argmax(distToLine)
    return idxOfBestPoint

def preprocess_rna_from_raw_counts(ad, 
                                   remove_mito_ribo = True,
                                   filter_genes_min_cells = 20,
                                   n_hvg = 1500, 
                                   n_pcs = 100,
                                   knee = None,
                                   neighbors = 30,
                                   umap_min_dist = 0.2):
    
    # remove mitochrondrial and ribosomal genes
    ad.var["mt"] = ad.var_names.str.startswith("MT-")
    # ribosomal genes
    ad.var["rb"] = ad.var_names.str.startswith(("RPS", "RPL"))
    if remove_mito_ribo:
        ad = ad[:,~((ad.var['mt']) | (ad.var['rb']))].copy()
        
    # filter genes
    sc.pp.filter_genes(ad, min_cells = filter_genes_min_cells)
    
    # normalize
    ad.layers['raw_counts'] = ad.X.copy()
    ad.layers['median'] = ad.layers['raw_counts'].copy()
    sc.pp.normalize_total(ad, layer='median')
    ad.layers['log'] = ad.layers['median'].copy()
    sc.pp.log1p(ad, layer='log')
    ad.X = ad.layers['log']
    
    # hvg
    sc.pp.highly_variable_genes(ad, n_top_genes=n_hvg)
    
    # PCA
    sc.tl.pca(ad, n_comps=n_pcs, use_highly_variable=True)
    ad.obsm['X_pca_max'] = ad.obsm['X_pca'].copy()
    if knee == None:
        curve = np.cumsum(ad.uns['pca']['variance_ratio'])
        knee = kneepoint(curve)
    ad.obsm['X_pca'] = ad.obsm['X_pca_max'][:, :knee]
    
    # UMAP
    sc.pp.neighbors(ad, method='umap', n_neighbors = neighbors, use_rep='X_pca', random_state = 5)
    sc.tl.umap(ad, min_dist = umap_min_dist, random_state=5)
    
    return ad

def coembed_anndata(sample_names, input_dir, outs_dir):

    # load in anndatas
    anndata_list = load_anndatas(sample_names, input_dir)
    
    comb_ad = combine_RNA_anndata(anndata_list)
    
    # preprocess
    comb_ad = preprocess_rna_from_raw_counts(comb_ad)
    
    # create out folder and write to it
    comb_ad.write(os.path.join(outs_dir, 'cr_merged_rna_ad.h5ad'))
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Co-embed RNA data from multiple filtered samples.')
    #parser.add_argument('--samples', '-s', nargs="*", type=str, help='Names of samples to include (should correspond to dir names in --data)')
    parser.add_argument('--data', '-i', type=str, help='Path to filtered sample output directory.')
    parser.add_argument('--out', '-o', type=str, help='Path to the output directory.')
    #parser.add_argument('--name', '-n', type=str, help='Name of this co-embedding.')
    args = parser.parse_args()
    input_dir = args.data
    outs_dir = args.out
    #name = args.name
    #samples = args.samples
    
    samples = [d.name for d in Path(input_dir).iterdir() if d.is_dir()]
    
    # You can also give samples by list
    #samples = ['56', '57', '58', '59', '60']
    
    samples.sort()
    print('processing samples: ' + str(samples))
    coembed_anndata(samples, input_dir, outs_dir)
