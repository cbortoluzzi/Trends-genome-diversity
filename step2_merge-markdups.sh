#!/bin/bash
#$ -N picard
#$ -cwd
#$ -pe thread 8
#$ -o output_%j.%N.txt
#$ -e error_%j.%N.txt
#$ -l h_rt=48:00:00
#$ -l s_rt=48:00:00
#$ -l h_vmem=40G
#$ -q long.q



bam1=$1		# E.g., BAZ-3086_TAATGCGC-ATAGAGGC-BH5TNKDSXX_L002.rg.sort.bam
bam2=$2		# E.g,, BAZ-3086_TAATGCGC-ATAGAGGC-BH5TNKDSXX_L003.rg.sort.bam


output=$(basename $bam1 | sed 's/_L002.rg.sort.bam//g')

mkdir -p results && mkdir -p results/bam

# Merge BAM files for a given sample
conda activate picard-2.26.2
picard MergeSamFiles I=$bam1 I=$bam2 VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true MERGE_SEQUENCE_DICTIONARIES=true O=results/bam/$output.rg.sort.merge.bam > results/bam/$output.stderr 2> results/bam/$output.stdout

# Mark duplicates
picard MarkDuplicates I=results/bam/$output.rg.sort.merge.bam O=$output.rg.sort.md.bam M=results/bam/$output.metrics CREATE_INDEX=true 2>> results/bam/$output.stderr >> results/bam/$output.stdout

# Flagstat analysis
conda activate samtools-1.9
samtools flagstat results/bam/$output.rg.sort.md.bam > results/bam/$output.flagstat 2>> results/bam/$output.stderr

# Check integrity of bam 
samtools quickcheck -v results/bam/*.rg.sort.md.bam >> bad_bams.fofn && echo 'all ok' || echo 'some files failed check, see bad_bams.fofn'

rm results/bam/$output.rg.sort.merge.bam
rm results/bam/$output.rg.sort.merge.bai

echo -e "Done!"
