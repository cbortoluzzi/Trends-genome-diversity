#!/bin/bash
#$ -N NJ
#$ -cwd
#$ -pe thread 8
#$ -o output_%j.%N.txt
#$ -e error_%j.%N.txt
#$ -l h_rt=48:00:00
#$ -l s_rt=48:00:00
#$ -l h_vmem=40G
#$ -q long.q


vcf=$1	# variants_gatk_SNP_hardFiltered_genFiltered_variants.1.39.vcf.gz 

# Construct a NJ tree from an identity-by-state distance relationship matrix
conda activate plink-1.90b6.21

# Convert VCF into PLINK format
plink --vcf $vcf --chr-set 39 --allow-extra-chr --recode --make-bed

# Obtain IBS matrix
plink --bfile plink --chr-set 39 --allow-extra-chr --maf 0.05 --geno 0.9 --distance

# Prepare input for phylip
perl 02.mdist2phylip.pl plink.mdist plink.fam inputPhylip

# Build NJ tree
conda activate phylip-3.697
phylip neighbor
