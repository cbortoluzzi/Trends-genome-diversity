#!/bin/bash
#$ -N liftover
#$ -cwd
#$ -pe thread 8
#$ -o output_%j.%N.txt
#$ -e error_%j.%N.txt
#$ -l h_rt=100:00:00
#$ -l s_rt=100:00:00
#$ -l h_vmem=60G
#$ -q long.q


target=$1	# This is the old genome in 2bit
query=$2	# This is the new genome in 2bit


# Run axtChain: chain together axt alignments
mkdir -p 2.axtChain/
for axt in 1.lastz/chr*/*.axt;do
        chain=$(basename $axt .axt)".chain"
        ./axtChain -linearGap=loose $axt $target $query 2.axtChain/$chain
done

# Run chainMergeSort: combine sorted files into larger sorted file
mkdir -p 3.chainMergeSort
./chainMergeSort 2.axtChain/*.chain > 3.chainMergeSort/all.chain
./chainSort 3.chainMergeSort/all.chain 3.chainMergeSort/all.sorted.chain

# Run chainPreNet - Remove chains that don't have a chance of being netted
mkdir -p 4.Net
./twoBitInfo $target $(basename $target .2bit)".sizes"
./twoBitInfo $query $(basename $query .2bit)".sizes"
./chainNet 3.chainMergeSort/all.sorted.chain $(basename $target .2bit)".sizes" $(basename $query .2bit)".sizes" 4.Net/all.net /dev/null
./netChainSubset 4.Net/all.net 3.chainMergeSort/all.chain galGal6ToGalGal7.chain

# Reformat chain
cat galGal6ToGalGal7.chain | grep -v '#' | gzip > galGal6ToGalGal7.over.chain.gz
rm galGal6ToGalGal7.chain

