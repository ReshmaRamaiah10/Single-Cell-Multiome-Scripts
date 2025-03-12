import scanpy as sc
import numpy as np
import pandas as pd
import anndata
from matplotlib import pyplot as plt
import seaborn as sns
import warnings
from anndata import AnnData
import sklearn
from sklearn.decomposition import PCA
from scipy.sparse import csr_matrix
import os

def subset_and_reprocess_rna(ad, 
                             cells_keep,
                             remove_mito_ribo = False,
                             filter_genes = False,
                             filter_genes_cutoff = 10,
                             n_hvg = 1500, 
                             n_pcs = 100,
                             knee = None,
                             neighbors = 30,
                             umap_min_dist = 0.3):

    # prep raw counts, obs and var
    try:
        count_mtx = ad[cells_keep,].layers['raw_counts']
    except KeyError: 
        count_mtx = ad[cells_keep,].X
    obs_df = ad[cells_keep,].obs
    var_df = ad[cells_keep,].var
    # join together
    sub_ad = anndata.AnnData(count_mtx, 
                          obs = obs_df, 
                          var = var_df)
    # remove mitochrondrial and ribosomal genes
    if remove_mito_ribo:
        sub_ad = sub_ad[:,~((sub_ad.var['mito']) | (sub_ad.var['ribo']))].copy()
    # remove lowly expressed genes
    if filter_genes:
        sc.pp.filter_genes(sub_ad, min_cells = filter_genes_cutoff)
    # normalize
    sub_ad.layers['raw_counts'] = sub_ad.X.copy()
    sub_ad.layers['median'] = sub_ad.layers['raw_counts'].copy()
    sc.pp.normalize_total(sub_ad, layer='median')
    sub_ad.layers['log'] = sub_ad.layers['median'].copy()
    sc.pp.log1p(sub_ad, layer='log')
    sub_ad.X = sub_ad.layers['log']
    # hvg
    sc.pp.highly_variable_genes(sub_ad, n_top_genes=n_hvg)
    # PCA
    sc.tl.pca(sub_ad, n_comps=n_pcs, use_highly_variable=True)
    sub_ad.obsm['X_pca_max'] = sub_ad.obsm['X_pca'].copy()
    if knee == None:
        curve = np.cumsum(sub_ad.uns['pca']['variance_ratio'])
        knee = kneepoint(curve)
    sub_ad.obsm['X_pca'] = sub_ad.obsm['X_pca_max'][:, :knee]
    # UMAP
    sc.pp.neighbors(sub_ad, 
                    method='umap', 
                    n_neighbors = neighbors, 
                    use_rep='X_pca', 
                    random_state = 5)
    sc.tl.umap(sub_ad, 
               min_dist = umap_min_dist, 
               random_state=5)
    # FDL
    sc.tl.draw_graph(sub_ad)
    return sub_ad

#computes kneepoint in PCA cumulative variance explained
def kneepoint(vec):
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

# Inputs
cr_ad = sc.read_h5ad('/athena/josefowiczlab/scratch/rer4011/projects/tori_atac_data/results/anndata/cr_merged_rna_ad.h5ad')
#metadata = pd.read_csv('/athena/josefowiczlab/scratch/rer4011/projects/tori_atac_data/results/metadata/01_merged_filtered_cells.csv',index_col=0)
metadata = pd.read_csv('/athena/josefowiczlab/scratch/rer4011/projects/tori_atac_data/results/metadata/01_merged_filtered_cells.csv')
out_ad = '/athena/josefowiczlab/scratch/rer4011/projects/tori_atac_data/results/anndata/manually_filtered1.ad'

metadata = metadata.astype(str)
metadata.index = metadata['orig.ident']+'#'+metadata['orig.barcode']
barcodes = metadata.index.tolist()

# Subset anndata
sub_ad = subset_and_reprocess_rna(cr_ad,barcodes)

# Merge the metadata
df_index_only = pd.DataFrame(index=sub_ad.obs.index)
merged_df = df_index_only.merge(metadata, left_index=True, right_index=True, how='left')
sub_ad.obs = merged_df
sub_ad.obs = sub_ad.obs.astype(str)

# Save the andata
sub_ad.write_h5ad(out_ad)
