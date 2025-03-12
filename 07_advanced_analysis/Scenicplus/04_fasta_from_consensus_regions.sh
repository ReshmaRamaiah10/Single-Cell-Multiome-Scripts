# Prepare fasta from consensus regions
# https://scenicplus.readthedocs.io/en/latest/human_cerebellum_ctx_db.html

REGION_BED="/data/niecr/cheongj/misc/results_seurat/scenicplus/output/consensus_peak_calling/consensus_regions.bed"
GENOME_FASTA="/home/cheongj/demultiplexing/hg38.fa"
CHROMSIZES="/data/niecr/cheongj/misc/qc_analysis_seurat/scenicplus/hg38.chrom.sizes"
DATABASE_PREFIX="MISC_pbmcpie_1kb_bg_with_mask"
SCRIPT_DIR="/data/niecr/cheongj/misc/qc_analysis_seurat/scenicplus/create_cisTarget_databases"

${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \
        ${GENOME_FASTA} \
        ${CHROMSIZES} \
        ${REGION_BED} \
        hg38.10x_brain.with_1kb_bg_padding.fa \
        1000 \
        yes
