import anndata
import os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import doubletdetection
import sys
import warnings

def subset_and_reprocess_rna(ad, 
                             cells_keep,
                             remove_mito_ribo = True,
                             filter_genes = True,
                             filter_genes_cutoff = 10,
                             n_hvg = 1500, 
                             n_pcs = 100,
                             knee = 45,
                             neighbors = 30,
                             umap_min_dist = 0.2):
    
    
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

def run_doubletdetection(ad, sample_col = 'sample', layer = None):
    ad.obs['doublet'] = np.nan
    ad.obs['doublet_score'] = np.nan
 
    # Calculate doublets on a per sample basis
    sample_names = ad.obs[sample_col].unique()
    for sample in sample_names:
        clf = doubletdetection.BoostClassifier()
 
        sub_ad = ad[ad.obs[sample_col] == sample, :]
        if layer == None:
            counts = sub_ad.X
        else: counts = sub_ad.layers[layer]
 
        warnings.filterwarnings('ignore')
        doublets = clf.fit(counts).predict(p_thresh=1e-7, voter_thresh=0.8)
        doublet_score = clf.doublet_score()
        warnings.filterwarnings('default')
 
        # Store doublets in adata
        ad.obs.loc[ad.obs[sample_col] == sample, 'doublet'] = doublets
        ad.obs.loc[ad.obs[sample_col] == sample, 'doublet_score'] = doublet_score
        
def calculate_qc_metrics(adata):
    # calculate qc metrics
    sc.pp.calculate_qc_metrics(adata, inplace = True)
    
    # ribo content
    RIBO_GENE_PATH = '/data/niecr/ramaiar/multiome_lucyrepo/analysis/mymodule/RB_genes_human'
    ribo_genes = pd.read_csv(RIBO_GENE_PATH, header = None)[0].tolist()
    ribo_genes = adata.var_names[adata.var_names.isin(ribo_genes)]
    adata.var['ribo'] = adata.var_names.isin(ribo_genes)
    row_sum_adata_ribo = np.sum(adata[:, ribo_genes].X.toarray(), axis=1)
    adata.obs['ribo_pct'] = row_sum_adata_ribo / adata.obs['total_counts'] * 100
    
    # mito content
    mt_genes = adata.var_names[adata.var_names.str.startswith('MT-')]
    adata.var['mito'] = adata.var_names.isin(mt_genes)
    row_sum_adata_mito = np.sum(adata[:, mt_genes].X.toarray(), axis=1)
    adata.obs['mito_pct'] = row_sum_adata_mito / adata.obs['total_counts'] * 100
    
    return adata
