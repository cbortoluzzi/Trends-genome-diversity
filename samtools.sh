#!/bin/bash
#$ -N samtools
#$ -cwd
#$ -pe thread 8
#$ -o output_%j.%N.txt
#$ -e error_%j.%N.txt
#$ -l h_rt=48:00:00
#$ -l s_rt=48:00:00
#$ -l h_vmem=40G
#$ -q long.q



bam=$1 # E.g., 194_ATTACTCG-CCTATCCT-BH5TNKDSXX.rg.sort.md.bam


# Calculate genome-wide coverage
conda activate samtools-1.9
samtools depth -a $bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > $bam.cov

# Calculate statistics
samtools stats --threads 32 $bam > $bam.stats


echo -e "Done!"
