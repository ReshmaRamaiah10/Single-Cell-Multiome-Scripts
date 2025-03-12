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

import doubletdetection
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
        
rna_ad = sc.read_h5ad('/athena/josefowiczlab/scratch/rer4011/projects/tori_atac_data/results/anndata/cr_merged_rna_ad.h5ad')
rna_ad

run_doubletdetection(rna_ad, sample_col = 'og_ident', layer = 'raw_counts')

rna_ad.write_h5ad('/athena/josefowiczlab/scratch/rer4011/projects/tori_atac_data/results/anndata/cr_merged_rna_ad.h5ad')
