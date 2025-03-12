#!/bin/bash

LOG_DIR='/data/niecr/cheongj/misc/qc_analysis_seurat/spectra/'

# processing script
PROC_SCRIPT='/data/niecr/cheongj/misc/qc_analysis_seurat/spectra/03_postprocess_spectra.py'

# specify run parameters
ADATA='/data/niecr/cheongj/misc/results_seurat/anndata_obj/03_spectra_annotated_rna_anndata.h5ad'
GS_DICT='/data/niecr/cheongj/misc/qc_analysis_seurat/spectra/misc_spectra_dict_resh_1.json'
MODEL_DIR='/data/niecr/cheongj/misc/results_seurat/spectra'
OBS_KEY='celltype_v1'

#### run postprocess script ####

H_DIRS=(${MODEL_DIR}/*)
# quick fix to run for only certain folders
#H_DIRS="/data/peer/sam/cytokine_central/models/liu21/spectra/v1/lambda001 /data/peer/sam/cytokine_central/models/liu21/spectra/v1/lambda005 /data/peer/sam/cytokine_central/models/liu21/spectra/v1/lambda01"

for i in ${H_DIRS[@]}
do
    MODEL_FILE=$(ls -1 "${i}"/*.pickle | xargs -n 1 basename)
    
    echo "running postprocess_spectra for $MODEL_FILE"

    # Submit the job
    bsub -J "spectra_post_${MODEL_FILE}" \
         -o "${LOG_DIR}/spectra_post_${MODEL_FILE}.%J.stdout" \
         -e "${LOG_DIR}/spectra_post_${MODEL_FILE}.%J.stderr" \
         -q gpuqueue \
         -n 2 \
         -R rusage[mem=100] \
         -W 24:00 \
         "source ~/.bashrc ; conda run -n multiome_winner_test3 python ${PROC_SCRIPT} ${ADATA} ${i}/${MODEL_FILE} NULL ${GS_DICT} ${i}/ ${OBS_KEY}"
done