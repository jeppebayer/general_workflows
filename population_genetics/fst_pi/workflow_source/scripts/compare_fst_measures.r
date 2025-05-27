#! /usr/bin/Rscript

### script for plotting fst results from popgen data

##############
## Packages
##############
# For plots
library(ggplot2)
# for colors
library(viridis)

check_names <- function(cols){
    # get all capitals
    cols <- gsub("aa", "A", cols)
    cols <- gsub("ae", "A", cols)
    cols <- gsub("oe", "O", cols)
    # check if there are semi-colons, if there is, not change should be done
    if( all(grepl(";", cols))){
        return(as.character(gsub(";", ":", cols)))    }
    # Split each part
    parts <- strcapture("^([a-zA-Z]{3})([^:]+):([a-zA-Z]{3})([^:]+)$", cols,
                    data.frame(id1A=character(), id1B=character(),
                               id2A=character(), id2B=character()))
    # Create short ID
    key <- paste(parts$id1A, parts$id2A, sep = ":")
    # Detect duplicates
    dup <- duplicated(key) | duplicated(key, fromLast = TRUE)
    # Check which side is ambiguous
    dup_id1A <- duplicated(parts$id1A) | duplicated(parts$id1A, fromLast = TRUE)
    dup_id2A <- duplicated(parts$id2A) | duplicated(parts$id2A, fromLast = TRUE)
    # Build output
    result <- ifelse(
        !dup, 
        paste(parts$id1A, parts$id2A, sep = ":"),
        ifelse(
            dup_id1A,
            paste(paste(parts$id1A, parts$id1B, sep = "."), parts$id2A, sep = ":"),
            paste(parts$id1A, paste(parts$id2A, parts$id2B, sep = "."), sep = ":")))
    result
}


##############
## Indputs
##############

fst_paths <- "/home/anneaa/EcoGenetics/general_workflow_outputs/population_genetics/fst/Collembola"
fst_files <- list.files(pattern = "mean", fst_paths, full.names = T, recursive = T)
fst_hud <- fst_files[grep("hudson", fst_files)]
fst_nei <- fst_files[-grep("hudson", fst_files)]; fst_nei <- fst_nei[-2]
fst_poolfstat <- list.files(pattern = "pairwise", fst_paths, full.names = T, recursive = T)
fst_poolfstat <- fst_poolfstat[grep("\\.fst", fst_poolfstat)]


##############
## plot Fst
##############

# fst
dat_hud <- lapply(fst_hud, function(x) read.table(x, header = T, sep = '\t', stringsAsFactors = F, check.names = FALSE))
dat_hud <- lapply(dat_hud, function(x) {rownames(x) <- x[,1]; x <- x[,-1]; x <- x["fst_mean",]})
names(dat_hud) <- fst_hud

dat_nei <- lapply(fst_nei, function(x) read.table(x, header = T, sep = '\t', stringsAsFactors = F, check.names = FALSE))
#dat_nei <- lapply(dat_nei, function(x) {colnames(x) <- x[2,]; x <- x[3,]})
names(dat_nei) <- fst_nei
#nam <- names(dat_nei[9]); dat_nei <- list(dat_nei[[9]]); names(dat_nei) <- nam

dat_poolf <- lapply(fst_poolfstat, function(x) read.table(x, header = T, sep = '\t', stringsAsFactors = F, check.names = FALSE))
dat_poolf <- lapply(dat_poolf, function(x) {xn <- rownames(x); x <- x[,"Fst Estimate"]; names(x) <- xn; return(as.data.frame(t(x)))})
names(dat_poolf) <- fst_poolfstat

fst_dat <- append(append(dat_hud, dat_nei), dat_poolf)

############
### get indv species
############

fst_dat_PogFla <- fst_dat[grep("Pogono", names(fst_dat))]
fst_dat_EntNic <- fst_dat[grep("Entomo", names(fst_dat))]
fst_dat_EntNic_GR <- fst_dat_EntNic[grep("grassland", names(fst_dat_EntNic))]
fst_dat_EntNic_CA <- fst_dat_EntNic[grep("conservation", names(fst_dat_EntNic))]

############
### get names similar
############

fst_dat_pog <- lapply(fst_dat_pog, function(x) { names(x) <- check_names(names(x)); x <- x[order(names(x))]; nm <- names(x); x <- t(x); rownames(x) <- nm; return(x)} )
str(fst_dat_pog)

fst_dat_EntNic_GR <- lapply(fst_dat_EntNic_GR, function(x) { names(x) <- check_names(names(x)); x <- x[order(names(x))]; nm <- names(x); x <- t(x); rownames(x) <- nm; return(x)} )
str(fst_dat_EntNic_GR)

fst_dat_EntNic_CA <- lapply(fst_dat_EntNic_CA, function(x) { names(x) <- check_names(names(x)); x <- x[order(names(x))]; nm <- names(x); x <- t(x); rownames(x) <- nm; return(x)} )
str(fst_dat_EntNic_CA)


############
### plot
############
plot_fun <- function(data_list){ 
    if (length(data_list) == 3){
        cols <- plasma(length(data_list))
        plot(data_list[[1]], col = cols[1], ylim = c(0,.5))
        for (dataset in seq(2,length(data_list))){
            #print(dataset)
            #print(data_list[[dataset]])
            points(data_list[[dataset]], col = cols[dataset], ylim = c(0,.5))
        }
    }
}

plot_fun(fst_dat_pog)
plot_fun(fst_dat_EntNic_GR)
plot_fun(fst_dat_EntNic_CA)
dev.off()

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


