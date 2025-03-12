#!/bin/bash

# List of file paths
file_paths=(
    "/athena/josefowiczlab/scratch/rer4011/projects/MAZ_andrew_data/raw_data/MAZ1"
    "/athena/josefowiczlab/scratch/rer4011/projects/MAZ_andrew_data/raw_data/MAZ2"
    "/athena/josefowiczlab/scratch/rer4011/projects/MAZ_andrew_data/raw_data/MAZ3"
    "/athena/josefowiczlab/scratch/rer4011/projects/MAZ_andrew_data/raw_data/MAZ4"
    "/athena/josefowiczlab/scratch/rer4011/projects/MAZ_andrew_data/raw_data/MAZ5"
    "/athena/josefowiczlab/scratch/rer4011/projects/MAZ_andrew_data/raw_data/MAZ6"
    "/athena/josefowiczlab/scratch/rer4011/projects/MAZ_andrew_data/raw_data/MAZ7"
    "/athena/josefowiczlab/scratch/rer4011/projects/MAZ_andrew_data/raw_data/MAZ8"
 )

for filepath in "${file_paths[@]}"; do
    samplename=${filepath##*/}
    sbatch --job-name="$samplename" 01_run_cellranger.sbatch "$filepath" "$samplename"  # Pass the sample name as the second argument
done
