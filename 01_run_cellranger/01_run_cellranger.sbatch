#!/bin/bash
#SBATCH -p scu-gpu
#SBATCH -c 8  # 8 cores
#SBATCH -t 0-144:00  # Running time of 144 hours
#SBATCH --mem 128G  # Memory request of 128 GB
#SBATCH -o /athena/josefowiczlab/scratch/rer4011/projects/MAZ_andrew_data/analysis/01_run_cellranger/%j
#SBATCH -e /athena/josefowiczlab/scratch/rer4011/projects/MAZ_andrew_data/analysis/01_run_cellranger/%j

# Add cellranger to PATH
export PATH=~/softwares/cellranger/cellranger-atac-2.1.0:$PATH

# Get filename from input
FILEPATH=${1}
samplename=${2}

# Print message indicating job submission
echo "Submitting job for sample: $samplename"

MOUSE_REF_PATH="/athena/josefowiczlab/scratch/rer4011/ref/refdata-cellranger-arc-mm10-2020-A-2.0.0"
WORKING_DIR="/athena/josefowiczlab/scratch/rer4011/projects/MAZ_andrew_data/cellranger_outs"

# Run cellranger atac
cd $WORKING_DIR

cellranger-atac count --id=$samplename \
                      --reference=$MOUSE_REF_PATH \
                      --fastqs=$FILEPATH \
                      --sample=$samplename \
                      --chemistry=ARC-v1
