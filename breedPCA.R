#!/usr/bin/env Rscript

library(gdsfmt)
library(argparse)
library(SNPRelate)
library(ggplot2)


set.seed(100)

args = commandArgs(trailingOnly=TRUE)

# The input VCF file can be bgzipped or not
vcf.fn <- args[1]
gds <- gsub(".vcf.gz$", ".gds", vcf.fn)

# Set color for plotting
#color <- '#998ec3' # Gasconne
color <- '#91cf60' # Barbezieux

# Reformat VCF 
#snpgdsVCF2GDS(vcf.fn, gds, method="biallelic.only")

# Open GDS file
genofile <- snpgdsOpen(gds)

# Add population information 
pop_code <- read.table(args[2], header=T, col.names=c("sample.id", "breed", "time")) # This tab delimited file should have information on: Sample code, breed name, and time point of sampling

# Principal component analysis - including LD filtering
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.5, autosome.only=TRUE, remove.monosnp=TRUE, missing.rate=10, maf=0.05)
snpset.id <- unlist(unname(snpset))
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2)

# Variance proportion (%)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

# Make a data.frame
samples <- data.frame("sample.id" = pca$sample.id)
dataframe_PCA <- data.frame(pop = merge(pop_code, samples, by="sample.id"), EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], EV3 = pca$eigenvect[,3], stringsAsFactors = FALSE)
output <- gsub(".vcf.gz$", "PCA.txt", vcf.fn)
write.table(dataframe_PCA, output, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# Plot PC1 vs PC2
dataframe_PCA$pop.time[dataframe_PCA$pop.time == 2015] <- 2013
plot <- ggplot(dataframe_PCA, aes(x=EV1, y=EV2, shape=as.factor(pop.time), size=2, alpha=0.2))+geom_point(colour=color)+theme_classic()+scale_color_brewer(palette="Accent")+xlab("PC1 (% of variance)")+ylab("PC2 (% of variance)")+theme(text = element_text(size = 20, family = 'sans'), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.text = element_text(color="black"))+theme(legend.position="none")
fig <- gsub(".vcf.gz$", ".pc1_vs_pc2.pdf", vcf.fn)
ggsave(fig, dpi=800, width=7, height=7)






