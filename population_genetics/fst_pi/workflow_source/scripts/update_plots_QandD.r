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
#args_are <- commandArgs(trailingOnly=T)

#fst_file <- args_are[1] 
    # fst_file <- "/home/anneaa/EcoGenetics/general_workflow_outputs/population_genetics/fst/Collembola/Entomobrya_nicoleti/EntNic_fst_allPos.fst"
#distance_file <- args_are[2]
    # distance_file = "/home/anneaa/EcoGenetics/general_workflows/population_genetics/fst/coordinates.tsv"
#outfile_ibd <- args_are[3]
#outfile_clad_neig <- args_are[4]
#outfile_clad_upgma <- args_are[5]
fst_paths <- "/home/anneaa/EcoGenetics/general_workflow_outputs/population_genetics/fst/Collembola"
fst_files <- list.files(pattern = "fst_allpairs_popoolation_mean.fst", fst_paths, full.names = T, recursive = T)
pi_file <- "/home/anneaa/EcoGenetics/general_workflow_outputs/population_genetics/pi/all_pi_estimates.txt"

##############
## plot Fst
##############

# fst
fst_dat <- lapply(fst_files, function(x) read.table(x, na.strings = "na", header = T, sep = '\t'))
names(fst_dat) <- gsub("conventional_agriculture", "K", gsub("conservation_agriculture", "CA", gsub("grassland", "GR", paste0(sapply(strsplit(sapply(strsplit(fst_files, "\\/"), "[[", 11), "_"), "[[", 1), "-", sapply(strsplit(sapply(strsplit(fst_files, "bola/"), "[[", 2), "\\/"), "[[", 2)))))
fst_dat <- lapply(fst_dat, function(x) x <- as.numeric(t(x[3,])))

str(fst_dat)

cols <- as.factor(sapply(strsplit(names(fst_dat),"-"), "[[", 1))
levels(cols) <- magma(length(levels(cols)), begin = .4, end = .98) 
png("/home/anneaa/EcoGenetics/people/anneaa/tests/intermediate_graphics_pi_fst/fst_allspecies.png", width = 7, height = 7, res = 600, units = "in")
par(oma = c(1,0,0,0))
boxplot(fst_dat, las = 2, col = as.character(cols), ylab = "Fst")
dev.off()

# cumbersome - for later
fst_dat <- lapply(fst_files, function(x) read.table(x, na.strings = "na", header = T, sep = '\t'))
names(fst_dat) <- gsub("conventional_agriculture", "K", gsub("conservation_agriculture", "CA", gsub("grassland", "GR", paste0(sapply(strsplit(sapply(strsplit(fst_files, "\\/"), "[[", 11), "_"), "[[", 1), "-", sapply(strsplit(sapply(strsplit(fst_files, "bola/"), "[[", 2), "\\/"), "[[", 2)))))
ignore_vec <- c("aeRoe", "FUR", "HaeJ", "HYF", "JHJ", "JYS", "KOS", "MYS", "RES", "SHJ", "K-DAS", "K-FLJ", "K-HVJ",
        "K-MoeJ", "K-ToeJ", "CA-BMJ", "CA-HVJ", "CA-HYS", "CA-JES", "CA-LES")
ignore_vec <- gsub("oe", "O", gsub("aa", "A", gsub("ae", "A", ignore_vec)))

# make function
remove_partial_match <- function(df, remove_vec){
        # Find columns to remove
    cols_to_remove <- sapply(df, function(col) {
        any(sapply(remove_vec, function(pattern) any(grepl(pattern, col))))
        })
        # Remove matched columns
    df_filtered <- df[ , !cols_to_remove]
}
# apply function across list
fst_dat$"IsoVir-CA"
lapply(fst_dat, function(x) ncol(x))
save_nm <- names(fst_dat)
fst_dat <- lapply(names(fst_dat), function(x) {
        values <- fst_dat[[x]]; 
        if(grepl("IsoVir", x)){ 
            values <- remove_partial_match(df = fst_dat[[x]], remove_vec = ignore_vec)}
        return(values)
        })
names(fst_dat) <- save_nm
lapply(fst_dat, function(x) ncol(x))
fst_dat$"IsoVir-CA"
# clean dfs
fst_dat <- lapply(fst_dat, function(x) x <- as.numeric(t(x[3,])))
str(fst_dat)

cols <- as.factor(sapply(strsplit(names(fst_dat),"-"), "[[", 1))
levels(cols) <- magma(length(levels(cols)), begin = .4, end = .98) 
png("/home/anneaa/EcoGenetics/people/anneaa/tests/intermediate_graphics_pi_fst/fst_IsoVirCleaned.png", width = 7, height = 7, res = 600, units = "in")
par(oma = c(1,0,0,0))
boxplot(fst_dat, las = 2, col = as.character(cols), ylab = "Fst")
dev.off()




# remove K field from isovir
fst_dat <- fst_dat[!grepl("IsoVir-K", names(fst_dat))]
cols <- as.factor(sapply(strsplit(names(fst_dat),"-"), "[[", 1))
levels(cols) <- magma(length(levels(cols)), begin = .4, end = .98) 
png("/home/anneaa/EcoGenetics/people/anneaa/tests/intermediate_graphics_pi_fst/fst_IsoVirCleaned_noK.png", width = 7, height = 7, res = 600, units = "in")
par(oma = c(1,0,0,0))
boxplot(fst_dat, las = 2, col = as.character(cols), ylab = "Fst")
dev.off()


# remove all IsoVir
fst_dat <- lapply(fst_files, function(x) read.table(x, na.strings = "na", header = T, sep = '\t'))
names(fst_dat) <- gsub("conventional_agriculture", "K", gsub("conservation_agriculture", "CA", gsub("grassland", "GR", paste0(sapply(strsplit(sapply(strsplit(fst_files, "\\/"), "[[", 11), "_"), "[[", 1), "-", sapply(strsplit(sapply(strsplit(fst_files, "bola/"), "[[", 2), "\\/"), "[[", 2)))))
fst_dat <- fst_dat[!grepl("IsoVir", names(fst_dat))]
fst_dat <- lapply(fst_dat, function(x) x <- as.numeric(t(x[3,])))
str(fst_dat)

png("/home/anneaa/EcoGenetics/people/anneaa/tests/intermediate_graphics_pi_fst/fst_noIsoVir.png", width = 7, height = 7, res = 600, units = "in")
cols <- as.factor(sapply(strsplit(names(fst_dat),"-"), "[[", 1))
levels(cols) <- magma(length(levels(cols)), begin = .4, end = .98) 
par(oma = c(1,0,0,0))
boxplot(fst_dat, las = 2, col = as.character(cols), ylab = "Fst")
dev.off()



# Plot pi
pi_dat <- read.table(pi_file, stringsAsFactors = T, header = T, sep = '\t', row.names = 1)
png("/home/anneaa/EcoGenetics/people/anneaa/tests/intermediate_graphics_pi_fst/Pi_allspecies.png", width = 7, height = 7, res = 600, units = "in")
cols <- as.factor(sapply(strsplit(colnames(pi_dat),"_"), "[[", 1))
levels(cols) <- magma(length(levels(cols)), begin = .4, end = .98) 
par(oma = c(1,0,0,0))
boxplot(pi_dat, las = 2, col = as.character(cols), ylab = "Pi")
dev.off()

# removing isovir data
pi_dat[which(rownames(pi_dat) %in% gsub("-","_", ignore_vec)), "IsoVir_CA"] <- NA
pi_dat[which(rownames(pi_dat) %in% gsub("-","_", ignore_vec)), "IsoVir_GR"] <- NA
pi_dat[which(rownames(pi_dat) %in% gsub("-","_", ignore_vec)), "IsoVir_K"] <- NA
png("/home/anneaa/EcoGenetics/people/anneaa/tests/intermediate_graphics_pi_fst/Pi_IsoVirCleaned.png", width = 7, height = 7, res = 600, units = "in")
cols <- as.factor(sapply(strsplit(colnames(pi_dat),"_"), "[[", 1))
levels(cols) <- magma(length(levels(cols)), begin = .4, end = .98) 
par(oma = c(1,0,0,0))
boxplot(pi_dat, las = 2, col = as.character(cols), ylab = "Pi")
dev.off()


