# Create cistarget databases
# https://scenicplus.readthedocs.io/en/latest/human_cerebellum_ctx_db.html

REGION_BED="/data/niecr/cheongj/misc/results_seurat/scenicplus/output/consensus_peak_calling/consensus_regions.bed"
GENOME_FASTA="/home/cheongj/demultiplexing/hg38.fa"
CHROMSIZES="/data/niecr/cheongj/misc/qc_analysis_seurat/scenicplus/hg38.chrom.sizes"
DATABASE_PREFIX="MISC_pbmcpie_1kb_bg_with_mask"
SCRIPT_DIR="/data/niecr/cheongj/misc/qc_analysis_seurat/scenicplus/create_cisTarget_databases"
OUT_DIR=""${PWD}""
CBDIR="${OUT_DIR}/aertslab_motif_colleciton/v10nr_clust_public/singletons"
FASTA_FILE="${OUT_DIR}/hg38.10x_brain.with_1kb_bg_padding.fa"
MOTIF_LIST="${OUT_DIR}/motifs.txt"

"${SCRIPT_DIR}/create_cistarget_motif_databases.py" \
    -f ${FASTA_FILE} \
    -M ${CBDIR} \
    -m ${MOTIF_LIST} \
    -o ${OUT_DIR}/${DATABASE_PREFIX} \
    --bgpadding 1000 \
    -t 20
