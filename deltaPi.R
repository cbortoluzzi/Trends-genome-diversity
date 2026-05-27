#!/usr/bin/env Rscript

library(ggplot2)


args = commandArgs(trailingOnly=TRUE)

# Read input file
df <- read.table(args[1], header=T)
df$Time_point[df$Time_point == 2015] <- 2013

# Plot genome-wide heterozygosity 
p <- ggplot(df, aes(x=Breed_name, y=Genome_heterozygosity, fill=as.factor(Time_point)))+geom_boxplot(position=position_dodge(0.8), notch=TRUE)+geom_dotplot(binaxis='y', stackdir='center',dotsize=0.5, position=position_dodge(0.8))+theme_classic()+scale_fill_brewer(palette="Accent")+ylim(0.003, 0.007)+theme(text = element_text(size = 20, family = 'sans'), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.text = element_text(color="black"))+theme(legend.position="none")+ylab("Heterozygosity (bp)")+xlab("")
ggsave('deltaPi.pdf', dpi=800, width=7, height=7)

# Calculate Wilcox test within breed but between time points
wilcox.test(subset(df$Genome_heterozygosity, df$Breed_name == "Barbezieux" & df$Time_point == 2003), subset(df$Genome_heterozygosity, df$Breed_name == "Barbezieux" & df$Time_point == 2013))
wilcox.test(subset(df$Genome_heterozygosity, df$Breed_name == "Gasconne" & df$Time_point == 2003), subset(df$Genome_heterozygosity, df$Breed_name == "Gasconne" & df$Time_point == 2013))

# Plot correlation between genome-wide heterozygosity and fraction genome in ROH
p <- ggplot(df, aes(x=Genome_heterozygosity, y=Genome_perc, colour=factor(Breed_name), shape=factor(Time_point), size=2, alpha=0.5))+geom_point()+theme_classic()+scale_color_brewer(palette="Accent")+theme(text = element_text(size = 20, family = 'sans'), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.text = element_text(color="black"))+theme(legend.position="none")+xlab("Heterozygosity (bp)")+ylab("Fraction of the genome covered by ROH")
ggsave('pi_vs_ROH.pdf', dpi=800, width=7, height=7)

# Calculate correlation between genome-wide heterozygosity and fraction genome in ROH
cor_test <- cor.test(df$Genome_heterozygosity, df$Genome_perc)

# Plot heterozygosity outside ROHs
p <- ggplot(df, aes(x=Breed_name, y=Heterozygosity_outsideROH, fill=as.factor(Time_point)))+geom_boxplot(position=position_dodge(0.8), notch=TRUE)+geom_dotplot(binaxis='y', stackdir='center',dotsize=0.5, position=position_dodge(0.8))+theme_classic()+scale_fill_brewer(palette="Accent")+ylim(0.003, 0.007)+theme(text = element_text(size = 20, family = 'sans'), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.text = element_text(color="black"))+theme(legend.position="none")+ylab("Heterozygosity outside ROH (bp)")+xlab("")
ggsave('pi_outsideROH.pdf', dpi=800, width=7, height=7)

# Calculate Wilcox test within breed but between time points
wilcox.test(subset(df$Heterozygosity_outsideROH, df$Breed_name == "Barbezieux" & df$Time_point == 2003), subset(df$Heterozygosity_outsideROH, df$Breed_name == "Barbezieux" & df$Time_point == 2013))
wilcox.test(subset(df$Heterozygosity_outsideROH, df$Breed_name == "Gasconne" & df$Time_point == 2003), subset(df$Heterozygosity_outsideROH, df$Breed_name == "Gasconne" & df$Time_point == 2013))

# Plot number of short, medium, and long ROH
p <- ggplot(df, aes(x=Breed_name, y=Num_short_ROH, fill=as.factor(Time_point)))+geom_boxplot(position=position_dodge(0.8))+geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(0.8), dotsize=0.5)+theme_classic()+scale_fill_brewer(palette="Accent")+theme(text = element_text(size = 20, family = 'sans'), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.text = element_text(color="black"))+theme(legend.position="none")+xlab("")+ylab("Number of short ROH")
ggsave('num_short_ROH.pdf', dpi=800, width=7, height=7)

p <- ggplot(df, aes(x=Breed_name, y=Num_medium_ROH, fill=as.factor(Time_point)))+geom_boxplot(position=position_dodge(0.8))+geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(0.8), dotsize=0.5)+theme_classic()+scale_fill_brewer(palette="Accent")+theme(text = element_text(size = 20, family = 'sans'), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.text = element_text(color="black"))+theme(legend.position="none")+xlab("")+ylab("Number of medium ROH")
ggsave('num_medium_ROH.pdf', dpi=800, width=7, height=7)

p <- ggplot(df, aes(x=Breed_name, y=Num_long_ROH, fill=as.factor(Time_point)))+geom_boxplot(position=position_dodge(0.8))+geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(0.8), dotsize=0.5)+theme_classic()+scale_fill_brewer(palette="Accent")+theme(text = element_text(size = 20, family = 'sans'), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.text = element_text(color="black"))+theme(legend.position="none")+xlab("")+ylab("Number of long ROH")
ggsave('num_long_ROH.pdf', dpi=800, width=7, height=7)

# Plot length of short, medium, long ROH (going to plot the sum in Mb)
p <- ggplot(df, aes(x=Breed_name, y=Sum_short_ROH/1000000, fill=as.factor(Time_point)))+geom_boxplot(position=position_dodge(0.8))+geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(0.8), dotsize=0.5)+theme_classic()+scale_fill_brewer(palette="Accent")+theme(text = element_text(size = 20, family = 'sans'), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.text = element_text(color="black"))+theme(legend.position="none")+xlab("")+ylab("Sum of short ROH (Mb)")
ggsave('sum_short_ROH.pdf', dpi=800, width=7, height=7)

p <- ggplot(df, aes(x=Breed_name, y=Sum_medium_ROH/1000000, fill=as.factor(Time_point)))+geom_boxplot(position=position_dodge(0.8))+geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(0.8), dotsize=0.5)+theme_classic()+scale_fill_brewer(palette="Accent")+theme(text = element_text(size = 20, family = 'sans'), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.text = element_text(color="black"))+theme(legend.position="none")+xlab("")+ylab("Sum of medium ROH (Mb)")
ggsave('sum_medium_ROH.pdf', dpi=800, width=7, height=7)

p <- ggplot(df, aes(x=Breed_name, y=Sum_long_ROH/1000000, fill=as.factor(Time_point)))+geom_boxplot(position=position_dodge(0.8))+geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(0.8), dotsize=0.5)+theme_classic()+scale_fill_brewer(palette="Accent")+theme(text = element_text(size = 20, family = 'sans'), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.text = element_text(color="black"))+theme(legend.position="none")+xlab("")+ylab("Sum of long ROH (Mb)")
ggsave('sum_long_ROH.pdf', dpi=800, width=7, height=7)

# Plot genomic or realised inbreeding
p <- ggplot(df, aes(x=Breed_name, y=FROH, fill=as.factor(Time_point)))+geom_boxplot(position=position_dodge(0.8))+geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(0.8), dotsize=0.5)+theme_classic()+scale_fill_brewer(palette="Accent")+theme(text = element_text(size = 20, family = 'sans'), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.text = element_text(color="black"))+theme(legend.position="none")+xlab("")+ylab("FROH")
ggsave('deltaFROH.pdf', dpi=800, width=7, height=7)

# Calculate Wilcox test within breed but between time points
wilcox.test(subset(df$FROH, df$Breed_name == "Barbezieux" & df$Time_point == 2003), subset(df$FROH, df$Breed_name == "Barbezieux" & df$Time_point == 2013))
wilcox.test(subset(df$FROH, df$Breed_name == "Gasconne" & df$Time_point == 2003), subset(df$FROH, df$Breed_name == "Gasconne" & df$Time_point == 2013))

