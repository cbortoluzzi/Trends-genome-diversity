#!/usr/bin/env python3
## Martijn Derks from Wageningen University, Netherlands and Maria Bernard from INRA, France
## Filter on average depth per sample set genotypes to missing

from collections import OrderedDict
import sys
import gzip
import argparse

class FILTER_VCF:

    def avg_depth_per_sample(self, table):
        self.sample_avg_depth_dic = {}
        depth_file = open(table,"r")
        for sample in depth_file:
            sample_id, sample_avg_depth = sample.strip().split()[0], sample.split()[1]
            self.sample_avg_depth_dic[sample_id] = float(sample_avg_depth.replace('X','')) ## Dictionary with avg depth per sample

    def filter_depth_per_sample(self, vcf, outfile):
        vcf_outfile = gzip.open(outfile,"wt")
        total_set_to_missing_un = 0
        total_set_to_missing_low = 0
        total_set_to_missing_high = 0
        nsites = 0
        format_idx = None
        depth_idx = None

        with gzip.open(vcf, 'rt', encoding = "ISO-8859-1") as vcf_reader:
            for record in vcf_reader:
                filtered_call = record.strip().split("\t")
                # variants line parsing
                if not record.startswith("#"):
                    record = record.strip().split("\t")
                    filtered_call = record[0:9]
                    sample_count = 0
                    nsites += 1
                    # identify depth in genotype format
                    if depth_idx is None:
                        depth_idx = record[format_idx].split(":").index("DP")
                    # parse samples genotype calling
                    for call in record[9:]:
                        # if known
                        if not call.startswith("."):
                            sample = samples[sample_count]
                            depth = call.split(":")[depth_idx]
                            avg_depth_times_2 = float(self.sample_avg_depth_dic[sample])*2.5
                            if depth == "." :
                                call = ".:.:."
                                total_set_to_missing_un += 1
                            elif int(depth) < 4 :
                                call = ".:.:."
                                total_set_to_missing_low += 1
                            elif float(depth) > avg_depth_times_2: ## Set to missing if depth > avg depth*2.5
                                call = ".:.:."
                                total_set_to_missing_high += 1
                        else:
                            pass #missing
                        filtered_call.append(call)
                        sample_count += 1
                # parsing genotype header line
                elif record.startswith("#CHROM"):
                    samples = record.split()[9:]
                    format_idx = record.split("\t").index("FORMAT")

                vcf_outfile.write("\t".join(filtered_call)+"\n")
        print( vcf, ": Avg set to missing per site: ", (total_set_to_missing_un + total_set_to_missing_low + total_set_to_missing_high)/nsites, "unknown depth : ", total_set_to_missing_un, "low depth : ", total_set_to_missing_low, "high depth : ", total_set_to_missing_high)

if __name__=="__main__":

    parser = argparse.ArgumentParser( description='Script to set genotype calls to missing if they do not contain coverage requirements')
    parser.add_argument("-v", "--vcf_file", help="VCF_file (compressed)")
    parser.add_argument("-b", "--bam_coverage", help="Tab delimited file with average depth per sample")
    parser.add_argument("-o", "--output_vcf", help="Output_VCF (gzipped)")
    args = parser.parse_args()
    F=FILTER_VCF()

    F.avg_depth_per_sample(args.bam_coverage)
    F.filter_depth_per_sample(args.vcf_file, args.output_vcf)


