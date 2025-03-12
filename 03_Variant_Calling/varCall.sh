#!/usr/bin/env bash

trim_galore --paired --basename $1 -j 8 ${1}*R{1,2}_001.fastq.gz
bwa-mem2 mem -t 8 hg38_bw2 ${1}*val* | \
  samtools sort -@ 8 - | \
  samtools view -@ 8 -q1 -b - > ${1}_sort.bam
picard MarkDuplicates VERBOSITY=WARNING \
  INPUT=${1}_sort.bam OUTPUT=${1}_rmDup.bam \
  REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT METRICS_FILE=tmp_${1}rmDup.log
picard DownsampleSam INPUT=${1}_rmDup.bam \
  OUTPUT=${1}_75pcReads.bam \
  RANDOM_SEED=50 PROBABILITY=0.75 VALIDATION_STRINGENCY=SILENT
samtools index ${1}_75pcReads.bam

picard DownsampleSam INPUT=${1}_rmDup.bam \
  OUTPUT=${1}_50pcReads.bam \
  RANDOM_SEED=50 PROBABILITY=0.50 VALIDATION_STRINGENCY=SILENT
samtools index ${1}_50pcReads.bam

picard DownsampleSam INPUT=${1}_rmDup.bam \
  OUTPUT=${1}_25pcReads.bam \
  RANDOM_SEED=50 PROBABILITY=0.25 VALIDATION_STRINGENCY=SILENT
samtools index ${1}_25pcReads.bam

bcftools mpileup --redo-BAQ --min-BQ 30 --per-sample-mF --annotate DP,AD -f /data/peer/lpaddock/data/ref/hg38/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \
  ${1}_rmDup.bam \
  ${1}_75pcReads.bam \
  ${1}_50pcReads.bam \
  ${1}_25pcReads.bam | \
  bcftools call --multiallelic-caller --variants-only -Ob > ${1}_bcftools.bcf

bcftools norm \
	-Ou \
	-m-any \
	 ${1}_bcftools.bcf \
	 | bcftools norm \
	 -Ov \
	 -f /data/peer/lpaddock/data/ref/hg38/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa > ${1}_bcftools.vcf

bcftools view -s ${1}_rmDup -c 1 ${1}_bcftools.vcf > ${1}_rmDup_variants.vcf
