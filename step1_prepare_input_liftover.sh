#!/bin/bash
#$ -N liftover
#$ -cwd
#$ -pe thread 8
#$ -o output_%j.%N.txt
#$ -e error_%j.%N.txt
#$ -l h_rt=100:00:00
#$ -l s_rt=100:00:00
#$ -l h_vmem=40G
#$ -q long.q


# Download old (softmasked) chicken reference genome - GRCg6a (Ensembl release 106)
wget http://ftp.ensembl.org/pub/release-106/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna_sm.toplevel.fa.gz

# Download new (softmasked) chicken reference genome - GRCg7b (Ensembl release 112)
wget https://ftp.ensembl.org/pub/release-112/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna_sm.toplevel.fa.gz

# Split old and new reference genome
mkdir -p chr_old && mkdir -p chr_new
# In the case of the new assembly, we will split it by creating 100 sequences
./faSplit sequence Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna_sm.toplevel.fa.gz 100 chr_new/ 
# In the case of the old assembly, we will split it by chromosome
./faSplit byname Gallus_gallus.GRCg6a.dna_sm.toplevel.fa.gz chr_old/

# Obtain 2bit file for both old and new reference genome
./faToTwoBit Gallus_gallus.GRCg6a.dna_sm.toplevel.fa.gz Gallus_gallus.GRCg6a.dna_sm.toplevel.2bit
./faToTwoBit Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna_sm.toplevel.fa.gz Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna_sm.toplevel.2bit
