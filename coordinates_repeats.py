#!/usr/bin/env python



# Author : Chiara Bortoluzzi



import argparse
from Bio import SeqIO
from itertools import groupby, count
from collections import defaultdict


parser = argparse.ArgumentParser(description = 'Obtain genomic coordinates of repeats')
parser.add_argument('--fasta', help = 'A soft-masked reference genome in FASTA format')



def coordinates_repeats(fasta):
	repeats = defaultdict(list)
	with open(fasta, "rt") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			chrom = record.id
			len_s = len(record.seq)
			for i in range(0, len_s):
				start = i
				nucleotide = record.seq[i]
				if nucleotide.islower():
					repeats[chrom].append(start)
	with open('coordinates_repeats.bed', 'w') as out:
		for key in repeats:
			coordinates = repeats[key]
			groups = groupby(coordinates, key=lambda item, c=count():item-next(c))
			tmp = [list(g) for k, g in groups]
			for i in tmp:
				first_pos = i[0]
				last_pos = i[-1] + 1
				out.write('{}\t{}\t{}\n'.format(key, first_pos, last_pos))



if __name__ == "__main__":
	args = parser.parse_args()
	repeats = coordinates_repeats(args.fasta)

