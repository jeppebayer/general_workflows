#! /usr/bin/Rscript

### script for plotting fst results from popgen data
#conda activate ecogen_neutral_diversity_wf; R
##############
## Packages
##############
# For plots
library(ggplot2)
# for colors
library(viridis)

check_names <- function(cols){
    # get all capitals
    #cols <- names(fst_dat_OrcCin_GR[[2]])
    # cols <- names(fst_dat_EntNic_GR[[3]])
    # cols <- names(fst_dat_PogFla[[3]])
    cols <- gsub("aa", "A", cols)
    cols <- gsub("ae", "A", cols)
    cols <- gsub("oe", "O", cols)
    cols <- gsub(".6J.", "J.", cols) # specific case
    cols <- gsub("-F-", "-", cols) # specific case
    cols <- gsub("\\.F\\.", "\\.", cols) # specific case
    # check if there are semi-colons, if there is, no change should be done
    if( any(grepl(";", cols))){
        cols <- gsub(";", ":", cols)
        if( any(grepl("\\.[0-9]{1}", cols)) ){# if the pattern of a . followed by a number is found
            cols <- gsub("\\.(?=[0-9]$)", "", cols, perl = TRUE)
            cols <- gsub("\\.(?=[0-9]:)", "", cols, perl = TRUE)
            return(cols)
        }
    }   
    if( any(grepl("-[0-9]{1}", cols)) ){# if the pattern of a - followed by a number is found
        cols <- gsub("-(?=[0-9]$)", "", cols, perl = TRUE)
        cols <- gsub("-(?=[0-9]:)", "", cols, perl = TRUE)
    }
        # Split each part
    parts <- strcapture("^([a-zA-Z]{3}).([^:]+):([a-zA-Z]{3}).([^:]+)$", cols,
                    data.frame(id1A=character(), id1B=character(),
                            id2A=character(), id2B=character()))
    if( any(is.na(parts))) {
        parts <- strcapture("^([a-zA-Z]{1}.[a-zA-Z]{3}).([^:]+):([a-zA-Z]{1}.[a-zA-Z]{3}).([^:]+)$", cols,
                    data.frame(id1A=character(), id1B=character(),
                            id2A=character(), id2B=character()))
    }
    if( any(is.na(parts))) {
        parts <- strcapture("^([a-zA-Z]{2}.[a-zA-Z]{3}).([^:]+):([a-zA-Z]{2}.[a-zA-Z]{3}).([^:]+)$", cols,
                    data.frame(id1A=character(), id1B=character(),
                                id2A=character(), id2B=character()))
    }  
    if( any(is.na(parts))) {
        parts <- strcapture("^([a-zA-Z]{2}.[a-zA-Z]{3}):([a-zA-Z]{2}.[a-zA-Z]{3})$", cols,
                    data.frame(id1A=character(), id2A=character()))
        parts$id1B <- NA
        parts$id2B <- NA
    }
    if( any(is.na(parts))) {
        parts <- strcapture("^([a-zA-Z]{3}):([a-zA-Z]{3})$", cols,
                    data.frame(id1A=character(), id2A=character()))
        parts$id1B <- NA
        parts$id2B <- NA
    }
    if( any(is.na(parts))) {
        parts <- strcapture( "^([a-zA-Z]{3})(?:\\.([A-Za-z0-9]+))?:([a-zA-Z]{3})(?:\\.([A-Za-z0-9]+))?$", 
            cols, data.frame(id1A=character(), id1B=character(), id2A=character(), id2B=character())
        )
    }
    if( any(is.na(parts))) {
        parts <- strcapture( "^([a-zA-Z]{3})(?:-([A-Za-z0-9]+))?:([a-zA-Z]{3})(?:-([A-Za-z0-9]+))?$", 
            cols, data.frame(id1A=character(), id1B=character(), id2A=character(), id2B=character())
        )
    }
    
    
    # Create short ID
    key <- paste(parts$id1A, parts$id2A, sep = ":")
    # Detect duplicates
    dup <- duplicated(key) | duplicated(key, fromLast = TRUE)

    # Subset to only duplicated rows
    dup_rows <- parts[dup, ]

    # Check ambiguity for id1A
    id1_ambiguous <- ave(dup_rows$id1B, dup_rows$id1A, FUN = function(x) length(unique(x)) > 1)

    # Check ambiguity for id2A
    id2_ambiguous <- ave(dup_rows$id2B, dup_rows$id2A, FUN = function(x) length(unique(x)) > 1)


    
    # Check which side is ambiguous
    #if(any(is.na(parts))){
    #    dup_id1Apart <- duplicated(parts[which(dup),]$id1A) | duplicated(parts[which(dup),]$id1A, fromLast = TRUE)
    #    dup_id2Apart <- duplicated(parts[which(dup),]$id2A) | duplicated(parts[which(dup),]$id2A, fromLast = TRUE)
    #    dup_id1A <- grepl("(A-Z){1}",parts$id1A)
    #    dup_id2A <- grepl("(A-Z){1}",parts$id2A)
    #    dup_id1A[which(dup)] <- dup_id1Apart
    #    dup_id2A[which(dup)] <- dup_id2Apart
    #}else{
    #    dup_id1A <- duplicated(parts$id1A) | duplicated(parts$id1A, fromLast = TRUE)
    #    dup_id2A <- duplicated(parts$id2A) | duplicated(parts$id2A, fromLast = TRUE)
    #}
    # Correct parts
    parts$id1A <- gsub("\\.", "-", parts$id1A)
    parts$id2A <- gsub("\\.", "-", parts$id2A)

    result <- character(length(cols))

    for (i in seq_along(cols)) {
        if (!dup[i]) {
            # No duplicates, just return compressed
            result[i] <- paste(parts$id1A[i], parts$id2A[i], sep = ":")
        } else {
            # In duplicated case, check ambiguity
            id1 <- parts$id1A[i]
            id2 <- parts$id2A[i]

            id1b_vals <- unique(parts$id1B[parts$id1A == id1])
            id2b_vals <- unique(parts$id2B[parts$id2A == id2])

            id1_amb <- length(id1b_vals) > 1
            id2_amb <- length(id2b_vals) > 1

            if (id1_amb && !id2_amb) {
            result[i] <- paste(paste(parts$id1A[i], parts$id1B[i], sep = "."), parts$id2A[i], sep = ":")
            } else if (!id1_amb && id2_amb) {
            result[i] <- paste(parts$id1A[i], paste(parts$id2A[i], parts$id2B[i], sep = "."), sep = ":")
            } else {
            # If both ambiguous or both not, include both
            result[i] <- paste(paste(parts$id1A[i], parts$id1B[i], sep = "."), paste(parts$id2A[i], parts$id2B[i], sep = "."), sep = ":")
            }
        }
    }


    # Build output
    #result <- ifelse(
    #    !dup, 
    #    paste(parts$id1A, parts$id2A, sep = ":"),
    #    ifelse(
    #        dup_id1A,
    #        paste(paste(parts$id1A, parts$id1B, sep = "."), parts$id2A, sep = ":"),
    #        paste(parts$id1A, paste(parts$id2A, parts$id2B, sep = "."), sep = ":")))
    result
}




plot_fun <- function(data_list){ 
    if (length(data_list) == 3){
        cols <- plasma(length(data_list))
        head_name_nr <- NROW(strsplit(names(data_list)[[1]], "/")[[1]])-2
        
        name_is <- sapply(strsplit(names(data_list)[[1]], "/"), "[[", head_name_nr)
        name_type <- sapply(strsplit(names(data_list)[[1]], "/"), "[[", head_name_nr+1)
        print(name_is)
        print(name_type)
        par(mar = c(8,4,1,1))
        plot(data_list[[1]], col = cols[1], ylim = c(0,1), xaxt = 'n', xlab = NA, ylab = "Fst")
        axis(1, labels = rownames(data_list[[1]]), at = seq(1, nrow(data_list[[1]])), las = 2, cex.axis = .5)
        title(main = gsub("_", " ", paste0(name_is,"\n", name_type)), line = -3)
        for (dataset in seq(2,length(data_list))){
            #print(dataset)
            #print(data_list[[dataset]])
            points(data_list[[dataset]], col = cols[dataset], ylim = c(0,1))
        }
        legend("topright", legend = c("Hudson", "Nei and Tajima", "Hivert - poolfstats"), col = cols, pch = 1)
    }
}

##############
## Indputs
##############

fst_paths <- "/home/anneaa/EcoGenetics/general_workflow_outputs/population_genetics/fst/Collembola"
fst_files <- list.files(pattern = "mean", fst_paths, full.names = T, recursive = T)
fst_hud <- fst_files[grep("hudson", fst_files)]
fst_nei <- fst_files[-grep("hudson", fst_files)]
fst_nei <- fst_nei[-grep("reg_mean", fst_nei)] ; fst_nei <- fst_nei[-2]
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
fst_dat_IsoVir <- fst_dat[grep("Isotoma", names(fst_dat))]
fst_dat_IsoVir_GR <- fst_dat_IsoVir[grep("grassland", names(fst_dat_IsoVir))]
fst_dat_IsoVir_CA <- fst_dat_IsoVir[grep("conservation", names(fst_dat_IsoVir))]
fst_dat_IsoVir_K <- fst_dat_IsoVir[grep("conven", names(fst_dat_IsoVir))]
fst_dat_OrcVir <- fst_dat[grep("Orchesella_vi", names(fst_dat))]
fst_dat_OrcVir_GR <- fst_dat_OrcVir[grep("grassland", names(fst_dat_OrcVir))]
fst_dat_OrcVir_CA <- fst_dat_OrcVir[grep("conservation", names(fst_dat_OrcVir))]
fst_dat_OrcCin <- fst_dat[grep("Orchesella_ci", names(fst_dat))]
fst_dat_OrcCin_GR <- fst_dat_OrcCin[grep("grassland", names(fst_dat_OrcCin))]

############
### get names similar
############


fst_dat_PogFla <- lapply(fst_dat_PogFla, function(x) { names(x) <- check_names(names(x)); x <- x[order(names(x))]; nm <- names(x); x <- t(x); rownames(x) <- nm; return(x)} )
str(fst_dat_PogFla)

fst_dat_EntNic_GR <- lapply(fst_dat_EntNic_GR, function(x) { names(x) <- check_names(names(x)); x <- x[order(names(x))]; nm <- names(x); x <- t(x); rownames(x) <- nm; return(x)} )
str(fst_dat_EntNic_GR)
fst_dat_EntNic_CA <- lapply(fst_dat_EntNic_CA, function(x) { names(x) <- check_names(names(x)); x <- x[order(names(x))]; nm <- names(x); x <- t(x); rownames(x) <- nm; return(x)} )
str(fst_dat_EntNic_CA)

fst_dat_IsoVir_CA <- lapply(fst_dat_IsoVir_CA, function(x) { names(x) <- check_names(names(x)); x <- x[order(names(x))]; nm <- names(x); x <- t(x); rownames(x) <- nm; return(x)} )
str(fst_dat_IsoVir_CA)
fst_dat_IsoVir_GR <- lapply(fst_dat_IsoVir_GR, function(x) { names(x) <- check_names(names(x)); x <- x[order(names(x))]; nm <- names(x); x <- t(x); rownames(x) <- nm; return(x)} )
str(fst_dat_IsoVir_GR)
fst_dat_IsoVir_K <- lapply(fst_dat_IsoVir_K, function(x) { names(x) <- check_names(names(x)); x <- x[order(names(x))]; nm <- names(x); x <- t(x); rownames(x) <- nm; return(x)} )
str(fst_dat_IsoVir_K)

fst_dat_OrcVir_GR <- lapply(fst_dat_OrcVir_GR, function(x) { names(x) <- check_names(names(x)); x <- x[order(names(x))]; nm <- names(x); x <- t(x); rownames(x) <- nm; return(x)} )
str(fst_dat_OrcVir_GR)
fst_dat_OrcVir_CA <- lapply(fst_dat_OrcVir_CA, function(x) { names(x) <- check_names(names(x)); x <- x[order(names(x))]; nm <- names(x); x <- t(x); rownames(x) <- nm; return(x)} )
str(fst_dat_OrcVir_CA)

fst_dat_OrcCin_GR <- lapply(fst_dat_OrcCin_GR, function(x) { names(x) <- check_names(names(x)); x <- x[order(names(x))]; nm <- names(x); x <- t(x); rownames(x) <- nm; return(x)} )
str(fst_dat_OrcCin_GR)



############
### plot
############


png("/home/anneaa/EcoGenetics/people/anneaa/tests/intermediate_graphics_pi_fst/fst_allspecies_pops.png", width = 7, height = 30, res = 600, units = "in")
par(mfrow = c(9,1))
plot_fun(fst_dat_EntNic_CA)
plot_fun(fst_dat_EntNic_GR)
plot_fun(fst_dat_IsoVir_CA)
plot_fun(fst_dat_IsoVir_GR)
plot_fun(fst_dat_IsoVir_K)
plot_fun(fst_dat_OrcVir_CA)
plot_fun(fst_dat_OrcVir_GR)
plot_fun(fst_dat_OrcCin_GR)
plot_fun(fst_dat_PogFla)
dev.off()










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


