#!/usr/bin/env python



# Author : Chiara Bortoluzzi



import argparse
from Bio import AlignIO
from collections import defaultdict


parser = argparse.ArgumentParser(description = 'Filter multiple sequence alignment by retaining alignment blocks where all species are present')
parser.add_argument('--maf', help = 'Multiple sequence alignment in MAF format containing a pairwise alignment between chicken and the reconstructed ancestral sequence')



def parse_MAF(maf_f):
	d = defaultdict(list)
	for multiple_alignment in AlignIO.parse(maf_f, "maf"):
		for seqrec in multiple_alignment:
			if seqrec.id.startswith('gallus_gallus.bGalGal1.mat.broiler.GRCg7b'):
				ref = seqrec.id
				refStart = seqrec.annotations["start"]
				refStrand = seqrec.annotations["strand"]
				refSize = seqrec.annotations["srcSize"]
				refLen = seqrec.annotations["size"]
				refSeq = seqrec.seq
				if refStrand == 1:
					newrefStrand = '+'
				else:
					newrefStrand = '-'
				d[ref, refStart, refLen, newrefStrand, refSize, refSeq]
			else:
				target = seqrec.id
				tStart = seqrec.annotations["start"]
				tStrand = seqrec.annotations["strand"]
				tSize = seqrec.annotations["srcSize"]
				tLen = seqrec.annotations["size"]
				tSeq = seqrec.seq
				if tStrand == 1:
					newtStrand = '+'
				else:
					newtStrand = '-'
				d[ref, refStart, refLen, newrefStrand, refSize, refSeq].append([target, tStart, tLen, newtStrand, tSize, tSeq])
	return d



def filter_MAF(d, output):
	with open(output, 'w') as out_bed:
		for key in d:
			for value in d[key]:
				for i in range(len(key[-1])):
					chrom = key[0].split('.')[-1]
					start = key[1] + i
					end = start + 1
					nuclref = key[5][i]
					nuclt = value[5][i]
					out_bed.write('{}\t{}\t{}\t{}\t{}\n'.format(chrom, start, end, nuclref, nuclt))



if __name__ == "__main__":
	args = parser.parse_args()
	output = args.maf.replace('.maf', '.bed')
	parser = parse_MAF(args.maf)
	filter = filter_MAF(parser, output)

