#!/bin/bash
#SBATCH -J cellranger_arc
#SBATCH -p scu-gpu
#SBATCH -c 8  # 8 cores
#SBATCH -t 0-144:00  # Running time of 144 hours
#SBATCH --mem 128G  # Memory request of 128 GB
#SBATCH -o /athena/josefowiczlab/scratch/rer4011/projects/MAZ_andrew_data/analysis/01_run_cellranger/%j
#SBATCH -e /athena/josefowiczlab/scratch/rer4011/projects/MAZ_andrew_data/analysis/01_run_cellranger/%j

# Add cellranger to PATH
export PATH=~/softwares/cellranger/cellranger-atac-2.1.0:$PATH

MOUSE_REF_PATH="/athena/josefowiczlab/scratch/rer4011/ref/refdata-cellranger-arc-mm10-2020-A-2.0.0"
ARC_CONFIG="/athena/josefowiczlab/scratch/rer4011/projects/MAZ_andrew_data/analysis/01_run_cellranger/02_cellranger_arc_config.csv"
WORKING_DIR="/athena/josefowiczlab/scratch/rer4011/projects/MAZ_andrew_data/cellranger_outs"

# Run cellranger atac
cd $WORKING_DIR

cellranger-atac aggr --id=combined \
                     --csv=$ARC_CONFIG \
                     --reference=$MOUSE_REF_PATH \
                     --normalize=none \
                     --dim-reduce=lsa