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



vcf=$1	# E.g., variants_gatk_all.vcf.gz


conda activate gatk4-4.2.4.0

# Generate a new VCF file containing the selected subset of variants
gatk SelectVariants -R ref/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa -V $vcf --select-type-to-include SNP -O results/gatk/variants_gatk_SNP.vcf.gz
gatk SelectVariants -R ref/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa -V $vcf --select-type-to-include INDEL -O results/gatk/variants_gatk_INDEL.vcf.gz

# Filter variants
gatk VariantFiltration -R ref/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa -V results/gatk/variants_gatk_SNP.vcf.gz --filter-expression "QD < 2.0; MQ < 40.0; FS > 60.0; SOR > 3.0; MQRankSum < -12.5; ReadPosRankSum < -8.0" --filter-name "SNP__HardFiltered" -O results/gatk/variants_gatk_SNP_hardfilter.vcf.gz
gatk VariantFiltration -R ref/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa -V results/gatk/variants_gatk_INDEL.vcf.gz --filter-expression "QD < 2.0; ReadPosRankSum < -20.0; InbreedingCoeff < -0.8; FS > 200.0; SOR > 10.0" --filter-name "INDEL_HardFiltered" -O results/gatk/variants_gatk_INDEL_hardfilter.vcf.gz

# Select variants
gatk SelectVariants -R ref/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa -V results/gatk/variants_gatk_SNP_hardfilter.vcf.gz --exclude-filtered -O results/gatk/variants_gatk_SNP_hardFiltered_variants.vcf.gz
gatk SelectVariants -R ref/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa -V results/gatk/variants_gatk_INDEL_hardfilter.vcf.gz --exclude-filtered -O results/gatk/variants_gatk_INDEL_hardFiltered_variants.vcf.gz

# Filter based on individual coverage
python3	filter_multi_sample_vcf.py --vcf_file results/gatk/variants_gatk_SNP_hardFiltered_variants.vcf.gz --bam_coverage sample_mean_cov.tsv --output_vcf results/gatk/variants_gatk_SNP_hardFiltered_genFiltered_variants.vcf.gz
python3 filter_multi_sample_vcf.py --vcf_file results/gatk/variants_gatk_INDEL_hardFiltered_variants.vcf.gz --bam_coverage sample_mean_cov.tsv --output_vcf results/gatk/variants_gatk_INDEL_hardFiltered_genFiltered_variants.vcf.gz

# Index VCF
conda activate htslib-1.14
tabix -p vcf results/gatk/variants_gatk_SNP_hardFiltered_genFiltered_variants.vcf.gz
tabix -p vcf results/gatk/variants_gatk_INDEL_hardFiltered_genFiltered_variants.vcf.gz


echo -e "Done!"
