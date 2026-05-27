#!/bin/bash
#$ -N gatk
#$ -cwd
#$ -pe thread 8
#$ -o output_%j.%N.txt
#$ -e error_%j.%N.txt
#$ -l h_rt=100:00:00
#$ -l s_rt=100:00:00
#$ -l h_vmem=40G
#$ -q long.q



gvcfs=`for f in results/gatk/*.g.vcf.gz;do echo "--variant $f";done | awk 'BEGIN { ORS = " " } { print }'`


# Merge several HaplotypeCaller GVCF files into a single GVCF with appropriate annotations
conda activate gatk4-4.2.4.0
gatk CombineGVCFs -R ref/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa $gvcfs -O results/gatk/cohort.g.vcf.gz 2> combineGVCFs.stderr > combineGVCFs.stdout

# Perform joint genotyping on one or more samples pre-called with HaplotypeCaller
gatk GenotypeGVCFs -R ref/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa -stand-call-conf 30 -V results/gatk/cohort.g.vcf.gz -O results/gatk/variants_gatk_all.vcf.gz 2>> combineGVCFs.stderr >> combineGVCFs.stdout

# Index VCF
conda activate htslib-1.14
tabix -p vcf results/gatk/variants_gatk_all.vcf.gz


echo -e "Done!"
