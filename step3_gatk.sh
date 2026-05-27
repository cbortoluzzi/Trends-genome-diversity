#!/bin/bash
#$ -N gatk
#$ -cwd
#$ -pe thread 8
#$ -o output_%j.%N.txt
#$ -e error_%j.%N.txt
#$ -l h_rt=120:00:00
#$ -l s_rt=120:00:00
#$ -l h_vmem=40G
#$ -q long.q



bam=$1	
vcf=$2	# Variants (including SNPs and indels) imported from dbSNP [remapped to GRCg7b]


output=$(basename $bam | sed 's/.rg.sort.md.bam//g')

mkdir -p results && mkdir -p results/gatk

# Generates recalibration table based on various user-specified covariates (such as read group, reported quality score, machine cycle, and nucleotide context).
conda activate gatk4-4.2.4.0
gatk BaseRecalibrator -R ref/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa -I $bam --known-sites $vcf -O results/gatk/$output.before.recal.table 2> results/gatk/$output.stderr > results/gatk/$output.stdout

# Apply a linear base quality recalibration model trained with the BaseRecalibrator tool.
gatk ApplyBQSR -bqsr results/gatk/$output.before.recal.table -I $bam -O results/gatk/$output.rg.sort.md.recal.bam -R ref/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa 2>> results/gatk/$output.stderr >> results/gatk/$output.stdout

# Call germline SNPs and indels via local re-assembly of haplotypes
gatk HaplotypeCaller -I results/gatk/$output.rg.sort.md.recal.bam -O results/gatk/$output.rg.sort.md.real.recal.g.vcf.gz -R ref/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa -ERC GVCF -mbq 10 --minimum-mapping-quality 30 2>> results/gatk/$output.stderr >> results/gatk/$output.stdout

# Index VCF
conda activate htslib-1.14
tabix -p vcf $output.rg.sort.md.real.recal.g.vcf.gz


echo -e "Done!"
