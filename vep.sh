#!/bin/bash
#$ -N VEP
#$ -cwd
#$ -pe thread 8
#$ -o output_%j.%N.txt
#$ -e error_%j.%N.txt
#$ -l h_rt=48:00:00
#$ -l s_rt=48:00:00
#$ -l h_vmem=40G
#$ -q long.q


# Obtain genomic coordinates of repeats from reference genome in FASTA format
python3 coordinates_repeats.py --fasta ../liftover/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna_sm.toplevel.fa

# Convert CADD score to BED format
zcat chCADD-scores/chCADD.tsv.gz | sed '1d' | awk 'OFS="\t"{print $1,$2-1,$2,$3,$4,$5}' > chCADD-scores/chCADD.bed

# Convert genomic coordaintes of CADD score to GalGal7
CrossMap bed ../liftover/galGal6ToGalGal7.over.chain.gz chCADD-scores/chCADD.bed chCADD-scores/chCADD.GalGal7.bed
rm -rf chCADD-scores/chCADD.bed

# Functional annotation of VCF file
for i in $(seq 1 39); do
	# Split CADD score by chromosome
	cat chCADD-scores/chCADD.GalGal7.bed | awk '{if($1 == "'$i'")print}' > chCADD-scores/chCADD.GalGal7.$i.bed

	# Run VEP to annotate bi-allelic SNPs
	conda activate ensembl-vep-110.1
	vep -i ../polarization/chromosome.$i.polarised.vcf --offline --species gallus_gallus --vcf --sift b -o chromosome.$i.polarised.vep.vcf --dir cache

	# Filter annotated VCF file based on call rate using VCFTools
	conda activate vcftools-0.1.16
	vcftools --vcf chromosome.$i.polarised.vep.vcf --max-missing 0.70 --recode --recode-INFO-all --out chromosome.$i.polarised.callrate.70.vep

	# Remove alleles that are found in repeats
	conda activate bedtools-2.30.0
	cat chromosome.$i.polarised.callrate.70.vep.recode.vcf | grep '#' > header.$i
	bedtools intersect -a chromosome.$i.polarised.callrate.70.vep.recode.vcf -b coordinates_repeats.bed -v > tmp.$i.vcf
	cat header.$i tmp.$i.vcf > chromosome.$i.polarised.callrate.70.no.repeats.vep.vcf
	rm -rf header.$i tmp.$i.vcf

	# Remove genes and transcripts that are not 1:1 orthologs between chicken and zebra finsh (genes downloaded from BioMart - Ensembl Release 112)
	python3 filter_orthologs.py --csv bioMart/chicken_zebrafinch_orthologs.csv --vcf chromosome.$i.polarised.callrate.70.no.repeats.vep.vcf

	# Assign the GERP score and chCADD score to all retained positions in the genome
	python3 assign_gerp_chCADD.py --bw GERP/gerp_conservation_scores.gallus_gallus.bGalGal1.mat.broiler.GRCg7b.bw --tsv chromosome.$i.polarised.callrate.70.no.repeats.orthologs.vep.tsv --CADD chCADD-scores/chCADD.GalGal7.$i.bed
done

rm -rf chromosome.*.polarised.callrate.70.vep.record.vcf
rm -rf chromosome.*.polarised.vep.vcf
rm -rf tmp.*.vcf
rm -rf header.*
