#!/usr/bin/env Rscript

#packages
library(poolfstat)


## inputs:
args_are <- commandArgs(trailingOnly = TRUE)

input_treemix_data_file <- args_are[1]
print(input_treemix_data_file)
# allele count data. populations in columns. the allele counts per population in the form of alleleA,alleleB. # nolint
# input_treemix_data_file <- "/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Beetles/Sitona_lineatus/intermediary_files/tmp/allele_freq/SitLin_Acount_allpops.count"
# input_treemix_data_file <- "/home/anneaa/EcoGenetics/people/anneaa/tests/allele_count_file.test"
# /home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Beetles/Sitona_lineatus/intermediary_files/tmp/allele_freq/SitLin_Acount_allpops.tmp1
# input_treemix_data_file <- "/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/grassland/Entomobrya_nicoleti/intermediary_files/tmp/allele_freq/EntNic_Acount_allpops.count"

position_file <- args_are[2]
print(position_file)
# positions of all alleles
# position_file <- "/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Beetles/Sitona_lineatus/intermediary_files/tmp/allele_freq/SitLin_# 
# position_file <- "/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Beetles/Sitona_lineatus/intermediary_files/tmp/allele_freq/SitLin_positions.tmp"
# position_file <- "/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/grassland/Entomobrya_nicoleti/intermediary_files/tmp/allele_freq/EntNic_positions.tmp"

#working_directory=$2
# working_directory=/home/anneaa/EcoGenetics/people/anneaa/tests
fst_output_dir <- args_are[3]
# fst_output_file=fst_file_all_pairs
# fst_output_dir <- "/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/grassland/Entomobrya_nicoleti/intermediary_files/tmp/fst"
# conda activate ecogen_neutral_diversity_wf
print(fst_output_dir)

species_short <- args_are[4]
#species_short <- "EntNic"

landscape_type <- args_are[5]
#landscape_type <- "grassland"


########################
### GET ALLELE counts  ###
########################

positions <- read.table(position_file, stringsAsFactors = F, header = F)

# positions <- positions[,c(1,3)]
colnames(positions) <- c("Chr", "Pos")
#allele_counts <- genotreemix2countdata(genotreemix.file = input_treemix_data_file, snp.pos = positions) # nolint: line_length_linter.
allele_counts <- genotreemix2countdata(genotreemix.file = input_treemix_data_file, snp.pos = positions, verbose = F) # nolint: line_length_linter.
	# filters can be applied
 
### modify pop names
prefix <- sub("\\..*", "", allele_counts@popnames)	 # Extract prefix before the dot
dup_prefix <- prefix[duplicated(prefix) | duplicated(prefix, fromLast = TRUE)] # Find duplicated prefixes
# Replace elements where prefix is not duplicated
allele_counts@popnames <- ifelse(prefix %in% dup_prefix, allele_counts@popnames, prefix)


allele_counts


########################
### Calculate overall FST  ###
########################

# using anova model by Weir & Cockerham 1984 and estimators by Weir 1996
fst_allelecounts <- computeFST(allele_counts)

#genome wide FST
fst_allelecounts$Fst
#> fst_allelecounts$Fst
# Entnic 41 pops:
#  Estimate bjack mean bjack s.e.    CI95inf    CI95sup 
# 0.1273034         NA         NA         NA         NA
write.table(fst_allelecounts$Fst, file = paste0(fst_output_dir, "/", species_short, "_", landscape_type, "_FST_poolfstat_overall_WeirCocker.tsv"), sep = "\t")
# Write out!

fst_allelecounts <- computeFST(allele_counts, nsnp.per.bjack.block = 1000)
#genome wide FST
fst_allelecounts$Fst
# Write out!
# entnic 41 pops
#   Estimate  bjack mean  bjack s.e.     CI95inf     CI95sup 
#0.127303391 0.127152007 0.002231927 0.122777431 0.131526584
write.table(fst_allelecounts$Fst, file = paste0(fst_output_dir, "/", species_short, "_", landscape_type, "_FST_poolfstat_blockJackn_WeirCocker.tsv"), sep = "\t")

# sliding window Fst
snps_in_window <- 50
fst_allelecounts <- computeFST(allele_counts, sliding.window.size = snps_in_window)

sliding_win_fst <- data.frame(fst_allelecounts$sliding.windows.fvalues$Chr, fst_allelecounts$sliding.windows.fvalues$CumMidPos/1e6, fst_allelecounts$sliding.windows.fvalues$MultiLocusFst)
colnames(sliding_win_fst) <- c("Chromosome", "Position_cummultive_Megabase", "Fst_multilocus")
head(sliding_win_fst)
# Write out!
write.table(sliding_win_fst, file = paste0(fst_output_dir, "/", species_short, "_", landscape_type, "_FST_overall_sliding_w", snps_in_window, "_poolfstat_WeirCocker.tsv"), sep = "\t")

# Write out!
dir.create(paste0(fst_output_dir, "/graphics"), showWarnings = F)
png(file = paste0(fst_output_dir, "/graphics/", species_short, "_", landscape_type, "_FST_poolfstat_overall_sliding_w", snps_in_window, "_WeirCocker_manhattan.png"), res = 600, units = "in", width = 12, height = 7)
plot(fst_allelecounts$sliding.windows.fvalues$CumMidPos/1e6, fst_allelecounts$sliding.windows.fvalues$MultiLocusFst,
	xlab = "Cumulated Position (in Mb)", ylab = paste0("Muli-locus Fst (", snps_in_window, ")"),
	col = as.numeric(as.factor(fst_allelecounts$sliding.windows.fvalues$Chr)), pch = 16, main = paste0(snps_in_window, " SNP windows"))
abline(h = fst_allelecounts$Fst, lty = 2)
dev.off()
# dashed line is overall fst, colors are alternating chromosomes

#without bundled
png(file = paste0(fst_output_dir, "/graphics/", species_short, "_", landscape_type, "_FST_poolfstat_overall_sliding_w", snps_in_window, "_WeirCocker_manhattan_bundRem.png"), res = 600, units = "in", width = 12, height = 7)
plot(fst_allelecounts$sliding.windows.fvalues[-which(fst_allelecounts$sliding.windows.fvalues$Chr == "bundled_sequences"),]$CumMidPos/1e6, fst_allelecounts$sliding.windows.fvalues[-which(fst_allelecounts$sliding.windows.fvalues$Chr == "bundled_sequences"),]$MultiLocusFst,
	xlab = "Cumulated Position (in Mb)", ylab = paste0("Muli-locus Fst (", snps_in_window, ")"),
	col = as.numeric(as.factor(fst_allelecounts$sliding.windows.fvalues[-which(fst_allelecounts$sliding.windows.fvalues$Chr == "bundled_sequences"),]$Chr)), 
	pch = 16, main = paste0(snps_in_window, " SNP windows"))
abline(h = fst_allelecounts$Fst, lty = 2)
dev.off()


###############################
### Calculate pairwise FST  ###
###############################

fst_allelecounts_pairwise <-compute.pairwiseFST(allele_counts, nsnp.per.bjack.block = 1000, verbose = F) #min.maf filters etc. can be employed
head(fst_allelecounts_pairwise@values)

# calc 95% ci: as jackknife mean +/- 1.96*S.E.
fst_allelecounts_pairwise@values$Fst_CI95l <- fst_allelecounts_pairwise@values$"Fst bjack mean" - 1.96 * fst_allelecounts_pairwise@values$"Fst bjack s.e."
fst_allelecounts_pairwise@values$Fst_CI95u <- fst_allelecounts_pairwise@values$"Fst bjack mean" + 1.96 * fst_allelecounts_pairwise@values$"Fst bjack s.e."

# print to file
head(fst_allelecounts_pairwise@values)
write.table(fst_allelecounts_pairwise@values, file = paste0(fst_output_dir, "/", species_short, "_", landscape_type, "_FST_pairwise_poolfstat", "_blockJackn_WeirCocker.fst"), sep = "\t")


# plot new type of fst plot. ordered with pop names
# heigth is dynamic
{
data <- data.frame(fst_allelecounts_pairwise@values$"Fst bjack mean", fst_allelecounts_pairwise@values$Fst_CI95l, fst_allelecounts_pairwise@values$Fst_CI95u, rownames(fst_allelecounts_pairwise@values))
colnames(data) <- c("Fst", "CI_l", "CI_u", "labels")
data <- data[order(data$Fst, decreasing = T),]
png(file = paste0(fst_output_dir, "/graphics/", species_short, "_", landscape_type, "_FST_pairwise_poolfstat", "_WeirCocker_fstplotT2.png"), res = 600, units = "in", width = 4, height = (0.03 * nrow(data) + 1))
par(mar = c(4, 8, 2, 2), cex = .5)  # space for y-axis labels
plot(
  data$Fst, 
  seq_along(data$Fst), 
  xlim = range(c(data$CI_l, data$CI_u)),  ylim = c(0,length(data$Fst)+1), #range(as.numeric(seq_along(data$Fst))),
  xlab = "Fst",   ylab = "",  pch = 19,   main = "Fst with 95% CIs", cex = .5, xaxs = "i", yaxs = "i", yaxt = 'n', col = "white"
)
# Horizontal error bars
segments(data$CI_l, seq_along(data$Fst), data$CI_u, seq_along(data$Fst), col = "gray50")
# points in black
points(data$Fst, seq_along(data$Fst), pch = 19, cex = .5)
# Y-axis labels
axis(2, at = seq_along(data$Fst), labels = data$labels, las = 2, cex.axis = .4, lwd.ticks = .5)
dev.off()
}


# regular fst plot (dot plot)
{
data <- data.frame(fst_allelecounts_pairwise@values$"Fst bjack mean", fst_allelecounts_pairwise@values$Fst_CI95l, fst_allelecounts_pairwise@values$Fst_CI95u, rownames(fst_allelecounts_pairwise@values))
colnames(data) <- c("Fst", "CI_l", "CI_u", "labels")
#data <- data[order(data$Fst, decreasing = T),]
png(file = paste0(fst_output_dir, "/graphics/", species_short, "_", landscape_type, "_FST_pairwise_poolfstat", "_WeirCocker_fstplotT1.png"), res = 600, units = "in", width = 11, height = 5)
par(mar = c(4.5, 4, 1, 1))  # space for y-axis labels
plot(data$Fst,
  ylab = "Fst",  pch = 19, col = "grey50", type = 'p', xlab = "Population pair", ylim = range(pretty(data$Fst))
)
dev.off()
}

# ordered but regular fst plot
{
data <- data.frame(fst_allelecounts_pairwise@values$"Fst bjack mean", fst_allelecounts_pairwise@values$Fst_CI95l, fst_allelecounts_pairwise@values$Fst_CI95u, rownames(fst_allelecounts_pairwise@values))
colnames(data) <- c("Fst", "CI_l", "CI_u", "labels")
data <- data[order(data$Fst, decreasing = T),]
png(file = paste0(fst_output_dir, "/graphics/", species_short, "_", landscape_type, "_FST_pairwise_poolfstat", "_WeirCocker_fstplotT1_ord.png"), res = 600, units = "in", width = 11, height = 5)
par(mar = c(4.5, 4, 1, 1))  # space for y-axis labels
plot(data$Fst,
  ylab = "Fst",  pch = 19, col = "grey50", type = 'p', xlab = "Population pair", ylim = range(pretty(data$Fst))
)
dev.off()
}


# plotting heatmap
{png(file = paste0(fst_output_dir, "/graphics/", species_short, "_", landscape_type, "_FST_pairwise_poolfstat", "_WeirCocker_heatmap.png"), res = 600, units = "in", width = 10, height = 10)
heatmap(fst_allelecounts_pairwise)
dev.off()}



# phylogeny
library(tidyr)
library(dplyr)
fst_vec <- fst_allelecounts_pairwise@values[,c(1)]
names(fst_vec) <- rownames(fst_allelecounts_pairwise@values)

fst_pairs <- data.frame(pair=rownames(fst_allelecounts_pairwise@values), fst=fst_allelecounts_pairwise@values$"Fst Estimate") %>%
  separate(pair, into = c("pop1", "pop2"), sep = ";")

# Get all unique populations
pops <- unique(c(fst_pairs$pop1, fst_pairs$pop2))

# Create empty matrix
fst_mat <- matrix(0, nrow=length(pops), ncol=length(pops),
                  dimnames = list(pops, pops))

# Fill the matrix
for (i in 1:nrow(fst_pairs)) {
  p1 <- fst_pairs$pop1[i]
  p2 <- fst_pairs$pop2[i]
  val <- fst_pairs$fst[i]
  fst_mat[p1, p2] <- val
  fst_mat[p2, p1] <- val
}

head(fst_mat)

# Convert to 'dist' object
fst_dist <- as.dist(fst_mat)
library(ape)
library("phangorn")
tree_upgma <- upgma(fst_dist) # does not exist
tree_nj <- nj(fst_dist)
{
png(file = paste0(fst_output_dir, "/graphics/", species_short, "_", landscape_type, "_FST_pairwise_poolfstat", "_WeirCocker_plylogeny_upgma.png"), res = 600, units = "in", width = 10, height = 10)
plot(tree_upgma, main="UPGMA Tree from FST")
dev.off()

png(file = paste0(fst_output_dir, "/graphics/", species_short, "_", landscape_type, "_FST_pairwise_poolfstat", "_WeirCocker_plylogeny_nj.png"), res = 600, units = "in", width = 10, height = 10)
plot(tree_nj, main="Neighbor Joining Tree from FST")
dev.off()
}
##############################
### Hierachical Fst
##############################

# groups <- 







