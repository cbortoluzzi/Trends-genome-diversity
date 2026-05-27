#!/bin/bash
#$ -N lastz
#$ -cwd
#$ -pe thread 8
#$ -o output_%j.%N.txt
#$ -e error_%j.%N.txt
#$ -l h_rt=100:00:00
#$ -l s_rt=100:00:00
#$ -l h_vmem=60G
#$ -q long.q


target=$1	# This is the old genome
query=$2	# This is the new genome

axt="chr"$(basename $target .fa)".vs."$(basename $query .fa)".axt"

# Run lastz
mkdir -p 1.lastz && mkdir -p 1.lastz/chr$(basename $target .fa)
echo -e $target vs $query
./lastz-1.04.00 --hspthresh=2200 --inner=2000 --ydrop=3400 --gappedthresh=10000 --scores=HoxD55 --chain --output=1.lastz/chr$(basename $target .fa)/$axt --format=axt --ambiguous=iupac ‑‑allocate:traceback=300M $target $query 2> $axt.stderr > $axt.stdout

