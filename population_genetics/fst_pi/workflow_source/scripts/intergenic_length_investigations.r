#! /usr/bin/Rscript

### script for plotting fst results from popgen data

##############
## Packages
##############
# For plots
library(ggplot2)
# for colors
library(viridis)

##############
## Indputs
##############
args_are <- commandArgs(trailingOnly=T)

intergenic_bed <- args_are[1] 
# intergenic_bed = "/home/anneaa/EcoGenetics/BACKUP/population_genetics/reference_genomes/collembola/Entomobrya_nicoleti/annotation/EG_EntNic_05092024_genomic.intergenic.bed"
working_directory <- args_are[2]
# working_directory = "/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/tmp"
##############
## Read in data
##############

# fst
intergenes <- read.table(intergenic_bed, stringsAsFactors = TRUE, na.strings = "na", header = F, sep = '\t')
head(intergenes)
tail(intergenes)
dim(intergenes)
colnames(intergenes) <- c("Scaffold", "Start", "End")

lengths <- intergenes[,3]-intergenes[,2]

# Where to output this?
pdf(paste0(working_directory, "/intergene_length_distribution.pdf"))
par(mfrow = c(2,1), oma = c(0,3,0,0), mar = c(4.7, 3.1, 1, 1))
plot(lengths, yaxt= 'n', xpd=NA, ylab = NA, main = NA, col = alpha("grey50", .6))
abline(h = mean(lengths), col = "red")
mtext(side=3, adj=.99, line=-1, paste("Mean =", round(mean(lengths))), col = "red", cex = .7)
mtext(side=3, adj=.99, line=-2, paste("Median =", round(median(lengths))), col = "blue", cex = .7)
abline(h = median(lengths), col = "blue")
title(ylab="Integenic lengths (bp)", mgp = c(4.5, 1, 0), xpd = NA)
axis(2, las = 2)
hist(lengths, breaks = 100, main = NA, xlab = "Integenic lengths (bp)", ylab = NA, yaxt = 'n')
axis(2, las = 2)
title(ylab="Frequency", mgp = c(4.5, 1, 0), xpd = NA)
dev.off()

# make new bed files with reduced intergenic regions based on percentage
head(intergenes)
for( percentage in c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)){
    lengths <- intergenes[,3]-intergenes[,2]
    remove_bases <- round( lengths * (percentage / 100) / 2)   # bases to remove on each end of intergene
    names <- c(colnames(intergenes)[1], paste0("Start_", percentage), paste0("End_", percentage))
    percentage_df <- data.frame(intergenes[,1], (intergenes$Start+remove_bases), intergenes$End-remove_bases)
    colnames(percentage_df) <- names
    print(head(percentage_df))
    write.table(percentage_df, file = paste0(working_directory, "/intergenic_", percentage,"_percent_rem.bed.temp"), sep = "\t")
}

# make new bed files with reduced intergenic regions based on number
head(intergenes)
for( number_remove in c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000)){
    lengths <- intergenes[,3]-intergenes[,2]
    #remove_bases <- round( lengths * (percentage / 100) / 2)   # bases to remove on each end of intergene
    names <- c(colnames(intergenes)[1], paste0("Start_", number_remove), paste0("End_", number_remove))
    percentage_df <- data.frame(intergenes[,1], (intergenes$Start+number_remove), intergenes$End-number_remove)
    colnames(percentage_df) <- names
    print(head(percentage_df))
    write.table(percentage_df, file = paste0(working_directory, "/intergenic_", number_remove,"_bases_rem.bed.temp"), sep = "\t")
}

## next use 
#pi_variant_pos_bed = /home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi/pi_allPops_variant_positions.bed
#al_freq_common_ps = /home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/allele_frequencies/Collembola/Entomobrya_nicoleti/allele_freq_allPops_all_positions.txt
## this is a common variant counter for all pops, instead of running on every single vcf?
#bedtools intersect -a pi_variant_pos_bed -b percentage.bed