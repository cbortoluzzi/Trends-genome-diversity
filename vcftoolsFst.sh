#!/bin/bash
#$ -N vcftools
#$ -cwd
#$ -pe thread 8
#$ -o output_%j.%N.txt
#$ -e error_%j.%N.txt
#$ -l h_rt=24:00:00
#$ -l s_rt=24:00:00
#$ -l h_vmem=40G
#$ -q long.q



vcf=$1	# variants_gatk_SNP_hardFiltered_genFiltered_variants.1.39.vcf.gz


# Calculate Fst in non-overlapping consecutive 50-kb windows - between breeds
vcftools --gzvcf $vcf --fst-window-size 50000 --weir-fst-pop gasconne.samples --weir-fst-pop barbezieux.samples --out barbezieux.vs.gasconne

# Calculate Fst in non-overlapping consecutive 50-kb windows - within breed, but between time points
vcftools --gzvcf $vcf --fst-window-size 50000 --weir-fst-pop gasconne.2003.samples --weir-fst-pop gasconne.2013.samples --out gasconne.2003.vs.2013
vcftools --gzvcf $vcf --fst-window-size 50000 --weir-fst-pop barbezieux.2003.samples --weir-fst-pop barbezieux.2013.samples --out barbezieux.2003.vs.2013
