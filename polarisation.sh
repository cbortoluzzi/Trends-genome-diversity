#!/bin/bash
#$ -N polarisation
#$ -cwd
#$ -pe thread 8
#$ -o output_%j.%N.txt
#$ -e error_%j.%N.txt
#$ -l h_rt=48:00:00
#$ -l s_rt=48:00:00
#$ -l h_vmem=40G
#$ -q long.q


chrom=$1

# This script assumes that the cactus alignment is already in path
# If it isn't then download it as: wget https://ftp.ensembl.org/pub/misc/compara/multi/hal_files/Fowl-10-way_20240131.hal

# Obtain a multiple sequence alignment in MAF format between the chicken and the ancestor of Galliformes (Anc26)
conda activate cactus-2.6.7
hal2maf --onlyOrthologs --refGenome gallus_gallus.bGalGal1.mat.broiler.GRCg7b --refSequence $chrom --targetGenomes Anc26 Fowl-10-way_20240131.hal chicken.vs.Anc26.$chrom.maf

# Transform MAF to BED format while filtering out alignment blocks where only the reference sequence is present
python3 filter_maf.py --maf chicken.vs.Anc26.$chrom.maf

# Split main VCF file into chromosome-level VCF files
conda activate bcftools-1.13
bcftools view -r $chrom ../gatk/variants_gatk_SNP_hardFiltered_genFiltered_variants.1.39.vcf.gz -Ov -o chromosome.$chrom.vcf --threads 20

# Polarise variants using information on ancestral sequence (only polarised bi-allelic SNPs will be retained)
python3 polarise_variants.py --bed chicken.vs.Anc26.$chrom.bed --vcf chromosome.$chrom.vcf

