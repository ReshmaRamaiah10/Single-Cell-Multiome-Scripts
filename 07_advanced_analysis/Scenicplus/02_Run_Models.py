# Run Models on Cistopic object
# Parallel LDA with MALLET
import os
import warnings
import pandas as pd
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk
warnings.simplefilter(action='ignore')
import pickle

# Project and output directories
projDir = '/data/niecr/cheongj/misc/results_seurat/scenicplus/'
out_dir = os.path.join(projDir, 'output/')

cistopic_obj = pickle.load(open(os.path.join(out_dir, "merged_cistopic_obj.pkl"), "rb"))

os.makedirs('scratch/ramaiar1/mallet/tutorial', exist_ok=True)

os.environ['MALLET_MEMORY'] = '200G'
from pycisTopic.lda_models import run_cgs_models_mallet
# Configure path Mallet
mallet_path="Mallet-202108/bin/mallet"
# Run models
models=run_cgs_models_mallet(
    cistopic_obj,
    n_topics= [30, 35, 40, 45, 50],
    n_cpu=12,
    n_iter=500,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    tmp_path="/scratch/ramaiar1/mallet/tutorial",
    save_path="/scratch/ramaiar1/mallet/tutorial",
    mallet_path=mallet_path,
)

pickle.dump(
    models,
    open(os.path.join(out_dir, "models.pkl"), "wb")
)
