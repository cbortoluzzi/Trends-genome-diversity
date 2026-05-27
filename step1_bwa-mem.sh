#!/bin/bash
#$ -N bwa-mem
#$ -cwd
#$ -pe thread 8
#$ -o output_%j.%N.txt
#$ -e error_%j.%N.txt
#$ -l h_rt=96:00:00
#$ -l s_rt=96:00:00
#$ -l h_vmem=60G
#$ -q long.q



ID=$1		# E.g., SAMEA115121650
read1=$2	# E.g., BAZ-3086_TAATGCGC-ATAGAGGC-BH5TNKDSXX_L002_R1.fastq.gz
read2=$3	# E.g., BAZ-3086_TAATGCGC-ATAGAGGC-BH5TNKDSXX_L002_R2.fastq.gz



# Prepare read groups for bwa-mem
header=$(gunzip -c $read1 | head -n 1)
sample=$(echo $(basename $read1 | cut -f1 -d'_'))
flowcell=$(echo $header | cut -f3 -d':')
lane=$(echo $header | head -n 1 | cut -f4 -d':')

output=$(echo $read1 | sed 's/_R1.fastq.gz//g' )


# Run bwa-mem2 mem
conda activate bwa-mem2-2.2.1
bwa-mem2 mem -M -t 32 -k 4 -R "@RG\tID:$flowcell.$lane\tPL:ILLUMINA\tSM:$sample\tPU:unknown\tCN:INRAE" ref/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa $read1 $read2 -o $output.sam 2> $output.stderr > $output.stdout

# Sort and convert SAM to BAM
conda activate samtools-1.9
samtools sort -@ 32 -T /tmp/ -o $output.rg.sort.bam -O BAM $output.sam 2>> $output.stderr >> $output.stdout
samtools index -@ 32 $output.rg.sort.bam 2>> $output.stderr >> $output.stdout

# Check integrity of BAM 
samtools quickcheck -v reads/*/*.rg.sort.bam > bad_bams.fofn && echo 'all ok' || echo 'some files failed check, see bad_bams.fofn'

rm $output.sam


echo -e "Done!"
