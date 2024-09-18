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

pattern <- args_are[1] 
# intergenic_bed = "/home/anneaa/EcoGenetics/BACKUP/population_genetics/reference_genomes/collembola/Entomobrya_nicoleti/annotation/EG_EntNic_05092024_genomic.intergenic.bed"
working_directory <- args_are[2]
# working_directory = "/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/tmp"
##############
## Read in data on variant in neutral regions
##############

files_variant_counts <- list.files(path = working_directory, pattern=paste0(".", pattern))

### Percentage removed from intergenic region:
files_variant_counts_perc <- files_variant_counts[grep("percent_rem", files_variant_counts)]

# read in data:
df <- data.frame()
for( file in files_variant_counts_perc ){
    # read in data
    data <- read.table(file, header = F)
    df <- data.frame(df, data)
    colnames(df) <- c(colnames(df, basename(file)))
}
head(df)

# Where to output this?
pdf(paste0(working_directory, "/snps_neutral_region_percentage_rem.pdf"))
plot(df[,1], yaxt= 'n', xpd=NA, ylab = NA, main = NA, col = alpha("grey50", .6), type = 'l', xlab = "Percentage removed")
axis(2, las = 2)
title(ylab="Variant count", mgp = c(4.5, 1, 0), xpd = NA)
par(new = TRUE)
for(each in length(2,ncol(df))){
    plot(df[,each], yaxt= 'n', xpd=NA, ylab = NA, main = NA, col = alpha("grey50", .6), type = 'l')
}
dev.off()


### N bases removed from intergenic region:
files_variant_counts_bases <- files_variant_counts[grep("_bases_rem", files_variant_counts)]

df <- data.frame()
for( file in files_variant_counts_bases ){
    # read in data
    data <- read.table(file, header = F)
    df <- data.frame(df, data)
    colnames(df) <- c(colnames(df, basename(file)))
}
head(df)

# Where to output this?
pdf(paste0(working_directory, "/snps_neutral_region_Nbases_rem.pdf"))
plot(df[,1], yaxt= 'n', xpd=NA, ylab = NA, main = NA, col = alpha("grey50", .6), type = 'l', xlab = "Bases removed")
axis(2, las = 2)
title(ylab="Variant count", mgp = c(4.5, 1, 0), xpd = NA)
par(new = TRUE)
for(each in length(2,ncol(df))){
    plot(df[,each], yaxt= 'n', xpd=NA, ylab = NA, main = NA, col = alpha("grey50", .6), type = 'l')
}
dev.off()




