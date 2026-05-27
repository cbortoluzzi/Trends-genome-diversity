#!/usr/bin/env python



# Author : Chiara Bortoluzzi


import gzip
import argparse


parser = argparse.ArgumentParser(description = 'Polarize variants based on ancestral and derived allele')
parser.add_argument('--bed', help = 'A BED file with information on reference and ancestral allele')
parser.add_argument('--vcf', help = 'A VCF file compressed with bgzip and indexed with tabix')



def parse_bed(bed_f):
	d = {}
	with open(bed_f) as f:
		for line in f:
			chrom, start, end, nuclchicken, nuclancestral = line.strip().split()
			chrom = int(chrom)
			start = int(start)
			nuclchicken = nuclchicken.upper()
			nuclancestral = nuclancestral.upper()
			d[chrom, start] = [nuclchicken, nuclancestral]
	return d


def parse_alignment(vcf_f, d):
	output_vcf = vcf_f.replace('.vcf', '.polarised.vcf')
	with open(output_vcf, 'w') as output:
		with open(vcf_f, 'r') as vcf_reader:
			for record in vcf_reader:
				if record.startswith('#'):
					header = record
					output.write('{}'.format(header))
				else:
					record = record.strip().split()
					chrom = int(record[0])
					pos = int(record[1])
					snpid = record[2]
					ref = record[3]
					alt = record[4]
					qual = record[5]
					filter = record[6]
					info = record[7]
					format = record[8]
					samples = '\t'.join(str(i) for i in record[9:])
					# We need to substract 1 to the position, because the VCF file is 1-based, whereas the BED file obtained from the MAF file is 0-based
					pos_zero_based = int(record[1]) - 1
					# We retain only bi-allelic
					if len(alt) == 1:
						try:
							reference = d[chrom, pos_zero_based][0]
							ancestral = d[chrom, pos_zero_based][1]
							# We want to double check that the reference in the MAF is the same as the reference in the VCF!
							assert reference == ref
							# Check if ancestral allele is equal to the reference or alternative allele
							if ancestral == ref:
								aa = "AA=" + str(ref)
								new_info = ';'.join(map(str, [aa, info]))
								output.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom, pos, snpid, ref, alt, qual, filter, new_info, format, samples))
							elif ancestral == alt:
								aa = "AA=" + str(alt)
								new_info = ';'.join(map(str, [aa, info]))
								output.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom, pos, snpid, ref, alt, qual, filter, new_info, format, samples))
						except KeyError:
							print (chrom, pos, "Ancestral not present")


if __name__ == "__main__":
	args = parser.parse_args()
	parse_ancestral = parse_bed(args.bed)
	pairwise_aln = parse_alignment(args.vcf, parse_ancestral)

