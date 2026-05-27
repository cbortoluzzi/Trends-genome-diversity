#!/usr/bin/env python



# Author : Chiara Bortoluzzi



import vcf
import argparse
import subprocess
from pathlib import Path



parser = argparse.ArgumentParser(description = 'Calculate genome-wide heterozygosity using a sliding window approach')
parser.add_argument('--vcf', help = 'A VCF file compressed with bgzip and indexed with tabix')
parser.add_argument('--bam', help = 'Aligned sequences in BAM format')
parser.add_argument('--tsv', help = 'Comma separated chromosome size file')
parser.add_argument('--minDP', help = 'Minimum depth to exclude a genotype [default: 4]', type = int, default = 4)
parser.add_argument('--maxDP', help = 'Maximum depth to exclude a genotype [default: 100]', type = int, default = 100)
parser.add_argument('--w', help = 'Window size [default = 10000 bp]', type = int, default = 10000)
parser.add_argument('--o', help = 'Output directory')



class Heterozygosity:

	def sequences_bam(self, tsv_f):
		"""
  		Example of a tsv_f:
    		CM046080.1,1,247284360
		...
		where CM046080.1 is the GenBank accession, 1 is the chromosome as integer, and 247284360 is the chromosome length
  		This function takes the GenBank accession and replaces it with an integer. This file should be prepared using
    		the genome assembly report file provided by NCBI.
  		"""
		self.mygenome = {}
		# Retain only autosomes and sex chromosomes
		with open(tsv_f) as f:
			for line in f:
				genbank_acc, chrom, length = line.strip().split(',')
				length = int(length)
				try:
					if isinstance(int(chrom), int):
						self.mygenome[genbank_acc] = [chrom, length]
				except ValueError:
					continue
		return self.mygenome


	def calculate_binned_heterozygosity(self, window, min_depth, max_depth, bam_f, vcf_f, output_file):
		for genbank_acc in self.mygenome:
			seq_length = self.mygenome[genbank_acc][1]
			for i in range(0, seq_length, window):
				start = i
				end = i + window
				self.bam_depth(genbank_acc, start, end, min_depth, max_depth, bam_f, vcf_f, window, output_file)


	def bam_depth(self, genbank_acc, start, end, min_depth, max_depth, bam_f, vcf_f, window, output_file):
		cov_sites = 0
		# Obtain the read depth of each site in the BAM file
		# Remember to load samtools before running this script ! 
		command = 'samtools depth -r %s:%d-%d %s' %(genbank_acc, start, end, bam_f)
		cmd = subprocess.check_output(command, shell = True).decode()
		outcmd = cmd.split('\n')
		for line in outcmd:
			if line:
				genbank_acc, position, depth = line.strip().split()
				# Filter sites based on read depth
				if int(depth) >= min_depth and int(depth) <= max_depth:
					cov_sites += 1
		self.heterozygosity(genbank_acc, start, end, cov_sites, vcf_f, window, output_file)


	def heterozygosity(self, genbank_acc, start, end, cov_sites, vcf_f, window, output_file):
		vcf_reader = vcf.Reader(filename=vcf_f, encoding = "ISO-8859-1")
		nhet = 0
		for record in vcf_reader.fetch(genbank_acc, start, end):
			for call in record.samples:
				nhet += record.num_het
		if cov_sites == window + 1:
			cov_sites = window
		try:
			SNPcount = round((window / cov_sites) * nhet, 3)
		except ZeroDivisionError:
			SNPcount = 0.0
		with open(output_file, 'a') as output_f:
			output_f.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(self.mygenome[genbank_acc][0], start, end, cov_sites, nhet, SNPcount))




if __name__ == "__main__":
	args = parser.parse_args()
	# create directory if it doesn't exist already
	path = Path(args.o)
	path.mkdir(parents=True, exist_ok=True)
	filename = Path(args.vcf).stem.replace('.vcf', '.het.' + str(args.w) + 'bp.txt')
	output_file = Path(path, filename)
	genome_wide_heterozygosity = Heterozygosity()
	genome_wide_heterozygosity.sequences_bam(args.tsv)
	genome_wide_heterozygosity.calculate_binned_heterozygosity(args.w, args.minDP, args.maxDP, args.bam, args.vcf, output_file)
