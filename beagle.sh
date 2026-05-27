#!/bin/bash
#$ -N beagle
#$ -cwd
#$ -pe thread 8
#$ -o output_%j.%N.txt
#$ -e error_%j.%N.txt
#$ -l h_rt=100:00:00
#$ -l s_rt=100:00:00
#$ -l h_vmem=60G
#$ -q long.q



vcf=$1


# Reformat input file (this is necessary because Beagle does not accept the way we set the missing genotypes)
conda activate htslib-1.14
zcat $vcf | perl -pe "s/\s\.:/\t.\/.:/g" | bgzip -c > input.beagle.vcf.gz
tabix -p vcf input.beagle.vcf.gz

for chr in $(seq 1 39);do
        # Phase haplotypes - necessary to call IBD
        java -jar beagle.08Feb22.fa4.jar gt=input.beagle.vcf.gz out=chromosome.$chr.phased chrom=$chr burnin=10 iterations=12 impute=true ne=100000 window=0.02 overlap=0.01 nthreads=16
        # Identify IBD segments
        java -jar refined-ibd.17Jan20.102.jar gt=chromosome.$chr.phased.vcf.gz out=chromosome.$chr.phased.IBD chrom=$chr length=0.02 window=0.06 lod=3 nthreads=12 trim=0.001
done





