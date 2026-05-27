#!/usr/bin/env Rscript

library(gdsfmt)
library(argparse)
library(SNPRelate)
library(ggplot2)

set.seed(1001)

args = commandArgs(trailingOnly=TRUE)

# The input VCF file can be bgzipped or not
vcf.fn <- args[1]

# Reformat VCF 
snpgdsVCF2GDS(vcf.fn, "allsamplesPCA.gds", method="biallelic.only")

# Open GDS file
genofile <- snpgdsOpen("allsamplesPCA.gds")

# Add population information 
pop_code <- read.table(args[2], header=T, col.names=c("sample.id", "breed", "time")) # This tab delimited file should have information on: Sample code, breed name, and time point of sampling

# Principal component analysis
pca <- snpgdsPCA(genofile, autosome.only=TRUE, remove.monosnp=TRUE, maf=0.05, missing.rate=10, num.thread=18)

# Variance proportion (%)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

# Make a data.frame
samples <- data.frame("sample.id" = pca$sample.id)
dataframe_PCA <- data.frame(pop = merge(pop_code, samples, by="sample.id"), EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], EV3 = pca$eigenvect[,3], stringsAsFactors = FALSE)
write.table(dataframe_PCA, "allsamplesPCA.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# Plot PC1 vs PC2
dataframe_PCA$pop.time[dataframe_PCA$pop.time == 2015] <- 2013
plot <- ggplot(dataframe_PCA, aes(x=EV1, y=EV2, colour=factor(pop.breed), shape=factor(pop.time), size=2, alpha=0.2))+geom_point()+theme_classic()+scale_color_brewer(palette="Accent")+xlab("PC1 (12.39% of variance)")+ylab("PC2 (4.76% of variance)")+theme(text = element_text(size = 20, family = 'sans'), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.text = element_text(color="black"))+theme(legend.position="none")
ggsave("allsamples_pc1_vs_pc2.pdf", dpi=800, width=7, height=7)

# Plot PC1 vs PC2
plot <- ggplot(dataframe_PCA, aes(x=EV1, y=EV3, colour=factor(pop.breed), shape=factor(pop.time), size=2, alpha=0.2))+geom_point()+theme_classic()+scale_color_brewer(palette="Accent")+xlab("PC1 (12.39% of variance)")+ylab("PC3 (3.41% of variance)")+theme(text = element_text(size = 20, family = 'sans'), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.text = element_text(color="black"))+theme(legend.position="none")
ggsave("allsamples_pc1_vs_pc3.pdf", dpi=800, width=7, height=7)




