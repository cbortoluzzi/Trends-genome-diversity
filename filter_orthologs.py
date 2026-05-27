#!/usr/bin/env python



# Author : Chiara Bortoluzzi


import pyBigWig
import argparse
from collections import defaultdict


parser = argparse.ArgumentParser(description = 'Remove variants whose gene/transcript is not 1:1 ortholog with zebra finch')
parser.add_argument('--csv', help = 'A BioMart file with information on chicken-zebra finch orthologs')
parser.add_argument('--vcf', help = 'A VCF file compressed with bgzip and indexed with tabix')




def one_to_one_ortholog(biomart):
	d = {}
	with open(biomart) as f:
		next(f)
		for line in f:
			line = line.strip().split(',')
			ortholog_type = line[-1]
			if ortholog_type == "ortholog_one2one":
				gene = line[0]
				transcript = line[2]
				d[gene, transcript] = ortholog_type
	return d


def filter_VCF(vcf_f, d):
	dict = defaultdict(lambda:defaultdict(list))
	with open(vcf_f, 'r') as vcf_reader:
		for record in vcf_reader:
			if record.startswith('#'):
				header = record
				line = header.strip()
				if line.startswith("#CHROM"):
					samples = line.split()[9:]
			else:
				record = record.strip().split()
				chrom = int(record[0])
				pos = int(record[1])
				ref = record[3]
				alt = record[4]
				info = record[7]
				format = record[8]
				allele_freq = info.split(';')[2].split('=')[1]
				genotypes = [item.split(':')[0] for item in  record[9:]]
				ancestral_allele = info.split(';')[0].replace('AA=','')
				annotation = info.split(';')[-1].split(',')
				for elem in annotation:
					# Format of VEP
					# Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SIFT
					elemt = elem.split('|')
					gene = elemt[4]
					transcript = elemt[6]
					try:
						ortholog = d[gene, transcript]
						if ortholog:
							geno_s = list(zip(samples, genotypes))
							for sample, geno in geno_s:
								if geno == '0/1' or geno == '0|1':
									dict[chrom, pos, ref, alt, ancestral_allele, allele_freq, elem]['het'].append(sample)
								else:
									# If the reference allele is the ancestral allele, then we choose 1/1
									if ref == ancestral_allele:
										if geno == '1/1' or geno == '1|1':
											dict[chrom, pos, ref, alt, ancestral_allele, allele_freq, elem]['hom'].append(sample)
									# Otherwise we choose 0/0
									elif alt == ancestral_allele:
										if geno == '0/0' or geno == '0|0':
											# And we recalculate the allele frequency
											AF = 1 - float(allele_freq)
											dict[chrom, pos, ref, alt, ancestral_allele, AF, elem]['hom'].append(sample)
					except KeyError:
						pass
	return dict




def filter_variants_using_RNA(dict, vcf_f):
#	bw_cerebellum = pyBigWig.open('RNAseq/bGalGal1.mat.broiler.GRCg7b.ENA.female_cerebellum.1.bam.bw')
	bw_heart = pyBigWig.open('RNAseq/bGalGal1.mat.broiler.GRCg7b.ENA.female_heart_muscle.1.bam.bw')
	bw_liver = pyBigWig.open('RNAseq/bGalGal1.mat.broiler.GRCg7b.ENA.female_liver.1.bam.bw')
	output = vcf_f.replace('.vep.vcf', '.orthologs.vep.tsv')
	with open(output, 'w') as out:
		for key in dict:
			chrom = str(key[0])
			# wiggle, bigWig, and bigBed files use 0-based half-open coordinates, which are also used by pyBigWig. So to access the value for the first base on chr1, one would specify the starting position as 0 and the end position as 1.
			# In our case, we need to substract one position to get the one we are interested in, because the VCF is 1-based.
			pos_zero_based = int(key[1]) - 1
			pos = int(key[1])
			#expr_cerebellum = bw_cerebellum.values(chrom, pos_zero_based, pos)[0]
			expr_heart = bw_heart.values(chrom, pos_zero_based, pos)[0]
			expr_liver = bw_liver.values(chrom, pos_zero_based, pos)[0]
			# Consider only positions that have an expression coverage of at least 200
			if expr_heart > 200 and expr_liver > 200:
				keys = '\t'.join(map(str, key))
				for geno in dict[key]:
					sample = ','.join(dict[key][geno])
					out.write('{}\t{}\t{}\n'.format(keys, geno, sample))




if __name__ == "__main__":
	args = parser.parse_args()
	one1one = one_to_one_ortholog(args.csv)
	filter_vcf = filter_VCF(args.vcf, one1one)
	expression = filter_variants_using_RNA(filter_vcf, args.vcf)

