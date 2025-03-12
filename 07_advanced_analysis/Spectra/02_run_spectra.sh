#!/bin/bash

# This script takes sevral days to finish running.
# Make sure to submit the job for 7 days

LOG_DIR='/data/niecr/cheongj/misc/qc_analysis_seurat/spectra/'

# processing script
PROC_SCRIPT='/data/niecr/cheongj/misc/qc_analysis_seurat/spectra/02_run_spectra.py'

# specify run parameters
ADATA='/data/niecr/cheongj/misc/results_seurat/anndata_obj/03_spectra_annotated_rna_anndata.h5ad'
GS_DICT='/data/niecr/cheongj/misc/qc_analysis_seurat/spectra/misc_spectra_dict_resh_1.json'
RESULT_DIR='/data/niecr/cheongj/misc/results_seurat/spectra'
LAMBDA_VALUES=( 0.1 0.05 0.01 0.0075 0.005 0.0025 0.001 )
OBS_KEY='celltype_v1'

for i in ${LAMBDA_VALUES[@]}
do
    echo "running spectra with lambda value $i"

    # Generate the bsub command
    bsub -J "spectra_${i}" \
         -o "${LOG_DIR}/spectra_${i}.%J.stdout" \
         -e "${LOG_DIR}/spectra_${i}.%J.stderr" \
         -q gpuqueue \
         -n 2 \
         -R rusage[mem=150] \
         -W 160:00 \
         "source ~/.bashrc ; conda run -n multiome_winner_test3 python ${PROC_SCRIPT} ${ADATA} ${GS_DICT} ${RESULT_DIR} ${i} ${OBS_KEY}"
done
