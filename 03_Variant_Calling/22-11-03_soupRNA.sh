#! /bin/bash -l

#SBATCH --partition=panda   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --job-name=soupRNA
#SBATCH --time=120:00:00   # HH/MM/SS
#SBATCH --mem=100G   # memory requested, units available: K,M,G,T
#SBATCH --output soupRNA-%j.out
#SBATCH --error soupRNA-%j.err

module load ~/anaconda3
source ~/.bashrc

conda activate souporcell

cd /athena/josefowiczlab/scratch/lup4006/IBD/souporcell

# ilc pre
singularity exec souporcell_latest.sif souporcell_pipeline.py \
    -i /athena/josefowiczlab/scratch/jic2016/ibd/ibd_ilc_pre/outs/gex_possorted_bam.bam \
    -b /athena/josefowiczlab/scratch/jic2016/ibd/ibd_ilc_pre/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
    -f /athena/josefowiczlab/scratch/jic2016/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
    -t 8 \
    -o ./22-11-03_ilc_pre_RNAsoup_out \
    -k 9
    
singularity exec Demuxafy.sif bash souporcell_summary.sh 22-11-03_ilc_pre_RNAsoup_out/clusters.tsv > 22-11-03_ilc_pre_RNAsoup_out/souporcell_summary.tsv

singularity exec Demuxafy.sif Assign_Indiv_by_Geno.R \
    -r /athena/josefowiczlab/scratch/lup4006/IBD/data/genotypes/22-10-28_IBD_genotypes_ilc_sorted.vcf \
    -c 22-11-03_ilc_pre_RNAsoup_out/cluster_genotypes.vcf \
    -o 22-11-03_ilc_pre_RNAsoup_out

# ilc post
singularity exec souporcell_latest.sif souporcell_pipeline.py \
-i /athena/josefowiczlab/scratch/jic2016/ibd/ibd_ilc_post/outs/gex_possorted_bam.bam \
-b /athena/josefowiczlab/scratch/jic2016/ibd/ibd_ilc_post/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
-f /athena/josefowiczlab/scratch/jic2016/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
-t 8 \
-o ./22-11-03_ilc_post_RNAsoup_out \
-k 9

singularity exec Demuxafy.sif bash souporcell_summary.sh 22-11-03_ilc_post_RNAsoup_out/clusters.tsv > 22-11-03_ilc_post_RNAsoup_out/souporcell_summary.tsv

singularity exec Demuxafy.sif Assign_Indiv_by_Geno.R \
    -r /athena/josefowiczlab/scratch/lup4006/IBD/data/genotypes/22-10-28_IBD_genotypes_ilc_sorted.vcf \
    -c 22-11-03_ilc_post_RNAsoup_out/cluster_genotypes.vcf \
    -o 22-11-03_ilc_post_RNAsoup_out

# pbmc pre
singularity exec souporcell_latest.sif souporcell_pipeline.py \
-i /athena/josefowiczlab/scratch/jic2016/ibd/ibd_pbmc_pre/outs/gex_possorted_bam.bam \
-b /athena/josefowiczlab/scratch/jic2016/ibd/ibd_pbmc_pre/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
-f /athena/josefowiczlab/scratch/jic2016/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
-t 8 \
-o ./22-11-03_pbmc_pre_RNAsoup_out \
-k 8

singularity exec Demuxafy.sif bash souporcell_summary.sh 22-11-03_pbmc_pre_RNAsoup_out/clusters.tsv > 22-11-03_pbmc_pre_RNAsoup_out/souporcell_summary.tsv

singularity exec Demuxafy.sif Assign_Indiv_by_Geno.R \
    -r /athena/josefowiczlab/scratch/lup4006/IBD/data/genotypes/22-10-28_IBD_genotypes_pbmc_sorted.vcf \
    -c 22-11-03_pbmc_pre_RNAsoup_out/cluster_genotypes.vcf \
    -o 22-11-03_pbmc_pre_RNAsoup_out

# pbmc post
singularity exec souporcell_latest.sif souporcell_pipeline.py \
-i /athena/josefowiczlab/scratch/jic2016/ibd/ibd_pbmc_post/outs/gex_possorted_bam.bam \
-b /athena/josefowiczlab/scratch/jic2016/ibd/ibd_pbmc_post/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
-f /athena/josefowiczlab/scratch/jic2016/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
-t 8 \
-o ./22-11-03_pbmc_post_RNAsoup_out \
-k 8

singularity exec Demuxafy.sif bash souporcell_summary.sh 22-11-03_pbmc_post_RNAsoup_out/clusters.tsv > 22-11-03_pbmc_post_RNAsoup_out/souporcell_summary.tsv

singularity exec Demuxafy.sif Assign_Indiv_by_Geno.R \
    -r /athena/josefowiczlab/scratch/lup4006/IBD/data/genotypes/22-10-28_IBD_genotypes_pbmc_sorted.vcf \
    -c 22-11-03_pbmc_post_RNAsoup_out/cluster_genotypes.vcf \
    -o 22-11-03_pbmc_post_RNAsoup_out

conda deactivate
