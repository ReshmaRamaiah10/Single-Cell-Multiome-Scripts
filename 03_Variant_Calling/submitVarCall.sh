#!/bin/bash -l

FASTQ_DIR="/data/niecr/cheongj/ibd/additional_genotype_fastq"

# get sample names
cd $FASTQ_DIR
samplenames=($(find . -type f -name "*R1*.fastq.gz" -exec basename {} \; | sed 's/_.*//'))

# set up job submission
for sample in "${samplenames[@]}"
do
    echo "Submitting job for sample: ${sample}"
    
    # Define job submission command here
    bsub -J vcf_${sample} \
         -o ${FASTQ_DIR}/vcf_${sample}.out \
         -e ${FASTQ_DIR}/vcf_${sample}.err \
         -q gpuqueue \
         -n 2 \
         -R "rusage[mem=100]" \
         -W 24:00 \
         -cwd ${FASTQ_DIR} \
         bash -l -c "source ~/.bashrc; conda activate varcaller; ./varCall.sh ${sample}"
done
