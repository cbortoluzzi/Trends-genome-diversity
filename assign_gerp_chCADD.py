#!/usr/bin/env python



# Author : Chiara Bortoluzzi



import pyBigWig
import argparse



parser = argparse.ArgumentParser(description = 'Assign the GERP and chCADD score to variants')
parser.add_argument('--bw', help = 'A GERP score file in bigWig format downloaded from Ensembl')
parser.add_argument('--CADD', help = 'A tab delimited file with chCADD information')
parser.add_argument('--tsv', help = 'A tab-delimited file with bi-allelic SNPs')



def parse_CADD(cadd_f):
	dict = {}
	with open(cadd_f) as f:
		for line in f:
			chrom, start, end, ref, alt, score = line.strip().split('\t')
			chrom = int(chrom)
			start = int(start)
			end = int(end)
			dict[chrom, start, end, ref, alt] = score
	return dict


def assign_gerp_and_CADD(bigwig, file, dict):
	output = file.replace('.vep.tsv', '.GERP.CADD.vep.tsv')
	bw_gerp = pyBigWig.open(bigwig)
	with open(output, 'w') as out:
		with open(file) as f:
			for line in f:
				line = line.strip().split()
				chrom = line[0]
				# wiggle, bigWig, and bigBed files use 0-based half-open coordinates, which are also used by pyBigWig. So to access the value for the first base on chr1, one would specify the starting position as 0 and the end position as 1.
				# In our case, we need to substract one position to get the one we are interested in, because the VCF is 1-based.
				pos_zero_based = int(line[1]) - 1
				pos = int(line[1])
				ref = line[2]
				alt = line[3]
				gerp_score = bw_gerp.values(chrom, pos_zero_based, pos)[0]
				info = '\t'.join(map(str, line))
				try:
					chrom = int(line[0])
					CADD_score = dict[chrom, pos_zero_based, pos, ref, alt]
					out.write('{}\t{}\t{}\n'.format(info, gerp_score, CADD_score))
				except KeyError:
					out.write('{}\t{}\t{}\n'.format(info, gerp_score, 'nan'))


if __name__ == "__main__":
	args = parser.parse_args()
	cadd = parse_CADD(args.CADD)
	gerp = assign_gerp_and_CADD(args.bw, args.tsv, cadd)

