"""
File: run_spectra_liu.py
Author(s): Suhani Balachandran, Sam Rose
Date: Oct 2023
Edited by: Reshma Ramaiah, July 31
Purpose: script to automate running Spectra on Josefowicz Lab MISC pilot
NOTES:
    * To print out command-line argument help menu: python3 run_spectra_lca.py -h
    * Results dir must have a subdirectory called lambda* where * is the lambda value without any decimal points. ie lambda = 0.1 --> lambda01
    * Runs with the default number of factors
    * Geneset must contain Ensembl IDs
"""

import pandas as pd
import numpy as np
from scipy import sparse
from copy import deepcopy
from collections import namedtuple
from datetime import date
import pickle
import scanpy as sc
import json
import anndata as ad
import sys
import os
import Spectra as spc
from Spectra import Spectra_util as spc_tl



# ---------------------------------------------------------------------
# STARTING MATERIALS
# ---------------------------------------------------------------------

# There should be processed anndata objects and gene set dictionaries created for the run
# these need to be specifified as well as an output directory, lambda value, and obs key for cell type annotations


# assign command line arguments to variables
adata_path = sys.argv[1]
geneset_path = sys.argv[2]
result_dir = sys.argv[3]
lambda_val = float(sys.argv[4])
obs_key = sys.argv[5]

if (sys.argv[1] == "-h" or sys.argv[1] == "--help"):
    print("The following arguments are required for run_spectra.py:\
          \n\t1. path to adata file\
          \n\t2. path to geneset dict file \
          \n\t3. path to results dir \
          \n\t4. lambda value\
          \n\t5. obs key for cell type annotations\
          \n\nExample usage: (assuming it is being called from within the /data/peer/suhani/glasner dir)\
          \npython run_spectra.py /data/adata.h5ad /data/peer/sam/cytokine-central/references/genesets/x.json /spectra/model 0.1\
          \n\nNOTES:\n\t* Results dir must have a subdirectory called lambda* where * is the lambda value without any decimal points. ie lambda = 0.1 --> lambda01\
          \n\t* Runs with the default number of factors\
          \n\t* Geneset must contain Ensembl IDs")
    

# ---------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------
## anndata
adata = sc.read_h5ad(adata_path)
# ensure dense matrix and log1p normalized counts
if sparse.issparse(adata.X):
    adata.X = adata.X.todense()

# gene set dictionary
with open(geneset_path, 'rb') as infile:
    annotations = json.load(infile)




# check gene set dictionary
## not necessary if the gene set dictionary is already in the correct format
#annotations = util.check_gene_set_dictionary(adata, annotations, obs_key = 'sublineage_coarse')

# set the number of global factors to 50 and the number of factors for each other key in annotations to their length + 3
# number will be just default of gene sets + 1
L_num = {}

for key in annotations.keys():
    if key == 'global':
        L_num[key] = 50
    else:
        L_num[key] = len(annotations[key]) + 3

# ---------------------------------------------------------------------
# RUN SPECTRA
# ---------------------------------------------------------------------

print("before spectra", flush=True)
# run spectra with specified lamda value
model = spc.est_spectra(adata = adata, gene_set_dictionary = annotations, 
                            cell_type_key = obs_key,
                            use_highly_variable = True, 
                            use_cell_types = True,
                            L = L_num,
                            use_weights = True, lam = lambda_val, 
                            delta=0.001,kappa = 0.00001, rho = 0.00001, 
                             n_top_vals = 50, 
                             label_factors = True,
                            num_epochs=5000)

# create output directory if it doesn't exist
out_dir = "{}/lambda{}".format(result_dir, str(lambda_val).replace(".", ""))

# Check if the directory exists
if not os.path.exists(out_dir):
    # If it doesn't exist, create it
    os.makedirs(out_dir)



# save the model
with open("{}/lambda{}_fdefault_{}.pickle".format(out_dir, str(lambda_val), str(date.today())), 'wb') as f:
    pickle.dump(model, f, pickle.HIGHEST_PROTOCOL)
    

print("done running", flush=True)



