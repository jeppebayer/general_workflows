#!/usr/bin/env Rscript

# conda activate migration_fsc
# cd /home/anneaa/EcoGenetics/general_workflows/population_genetics/fst_pi/configurations/collembola/Entomobrya_nicoleti

library(scales) # for color alpha function
library(tidyr)  # for data rearrangements
library(viridis)    # for colorblind friendly colors
library(magick) # for reading images
library(grid)   # for plotting image as raster
library(visreg)
#args_are <- commandArgs(trailingOnly = TRUE)


#path_to_search = args_are[1]
path_to_search = "/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/2dSFS/"
#pattern_is = args_are[2]
pattern_is = "bestlikelyhoods"
#name_append = args_are[3]
name_append = "EntNic_all_dat"

#list_of_files <- list.files(pattern="bestlhoods" , path="/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/2dSFS/2dSFS_EntNic_aaRJ_C225_vs_EntNic_aeRoe-C36_1_vs_2/", recursive = TRUE, full.names = TRUE)
# workaround before inputs are ready, and if input is unable to be collected and passed (due to length of argument)
list_of_files <- list.files(pattern = pattern_is , path = path_to_search, recursive = TRUE, full.names = TRUE)
NROW(list_of_files)
# remove pops!
# double check, they should already be removed in the pipeline. Better safe than sorry.
list_of_files1 <- list_of_files[grep("JEJ-C119|KoeJ-C212|VAJ-C24|HHJ-C162|SBJ-C52|LVJ-C76|ULJ-C122|SKJ-C16", list_of_files, invert = T)]
NROW(list_of_files1)
# list_of_files[grep("JEJ-C119|KoeJ-C212|VAJ-C24|HHJ-C162|SBJ-C52|LVJ-C76|ULJ-C122|SKJ-C16", list_of_files)]
#list_of_files1[grep("aaRJ",list_of_files1)][grep("600", list_of_files1[grep("aaRJ",list_of_files1)])]


    # only files from one pop comparison
#file <- "/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/2dSFS/2dSFS_EntNic_aaRJ_C225_vs_EntNic_aeRoe-C36_1_vs_2//900/2dSFS_EntNic_aaRJ_C225_vs_EntNic_aeRoe-C36_1_vs_2_900migDiv_sfs2d/2dSFS_EntNic_aaRJ_C225_vs_EntNic_aeRoe-C36_1_vs_2_900migDiv_sfs2d.bestlhoods"
df <- data.frame()
for(file in list_of_files1){
    datas <- read.table(file, header = T)
    name <- basename(file)
    count <- NROW(unlist(strsplit(name, "_")))
    migration_divide <- unlist(strsplit(name, "_"))[count-1]
    migration_divide <- as.integer(sapply(strsplit(migration_divide, "mig"), "[[", 1))
    #migration_divide <- migration_divide/3
    # get pop pair
    pop_pair <- gsub("2dSFS_", "", name)
    pop_pair <- gsub("_[0-9]*_vs_[0-9]*_[0-9]*migDiv_sfs2d.bestlikelyhoods", "", pop_pair)
    
    # make wide format
    datas1 <- data.frame(pop_pair, migration_divide, datas)
    df <- rbind(df, datas1)
}

str(df)
df_sort <- df[order(df[,1], df[,2], method = "auto"),]
nrow(df_sort)
#[1] 40210
df_sort_alldat <- df_sort
# remove pairs with split time below 7000
#df_sort <- df_sort[which(df_sort$TDIV > 7000),]
dat_rem = 14000
df_sort <- df_sort[which(df_sort$migration_divide <= dat_rem),]
dat_nrow <- nrow(df_sort); print(dat_nrow)
# 34744
#[1] 27143 (when 14000 generations)


df_sort <- df_sort[which(df_sort$TDIV > dat_rem),]
nrow(df_sort)
#[1] 27143 (when 14000 generations)
nrow(df_sort)/dat_nrow
#1] 0.7230313
# [1] 0.6750311 ( when 14000 genestaioN)

df_sort$migration_divide <- df_sort$migration_divide/3

# migration curve median

plotot <- function(){
    mean_frame <- aggregate(df_sort, list(df_sort$migration_divide), FUN = mean, na.rm = T)
    median_frame <- aggregate(df_sort, list(df_sort$migration_divide), FUN = median, rm.na = T)
    plot(mean_frame$Group.1, mean_frame$mig_rec, type = "p", ylim = c(min(c(mean_frame$mig_rec, mean_frame$mig_anc)), .0001), 
        ylab = "Gene flow rate", xlab = "Years bp")
    points(mean_frame$Group.1, mean_frame$mig_anc, col = "red", ylim = c(min(c(mean_frame$mig_rec, mean_frame$mig_anc)), .0001))
    points(median_frame$Group.1, median_frame$mig_rec, pch = 2, col = "black", ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), .0001))
    points(median_frame$Group.1, median_frame$mig_anc, pch = 2, col = "red", ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), .0001))
    legend("topright", legend=c("Ancestral migration", "Recent migration"), col = c("red", "black"), pch = 1, bty = "n")
    legend("topright", legend=c("Mean", "Median"), col = "black", pch = c(1,2), bty = "n", inset = c(.15,.09))
}
png(file = paste0("/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/median_migration_curve_", name_append,"_",dat_rem, ".png"), res = 800, units = "in", height = 7, width = 7)
plotot()
dev.off()

# migration curve median zoomed
plotot <- function(){
    mean_frame <- aggregate(df_sort, list(df_sort$migration_divide), FUN = mean, rm.na = T)
    median_frame <- aggregate(df_sort, list(df_sort$migration_divide), FUN = median, rm.na = T)
    plot(mean_frame$Group.1, mean_frame$mig_rec, type = "p", ylim = c(min(c(mean_frame$mig_rec, mean_frame$mig_anc)), .00002), 
        ylab = "Gene flow rate", xlab = "Years bp")
    points(mean_frame$Group.1, mean_frame$mig_anc, col = "red", ylim = c(min(c(mean_frame$mig_rec, mean_frame$mig_anc)), .00002))
    points(median_frame$Group.1, median_frame$mig_rec, pch = 2, col = "black", ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), .00002))
    points(median_frame$Group.1, median_frame$mig_anc, pch = 2, col = "red", ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), .00002))
    legend("topright", legend=c("Ancestral migration", "Recent migration"), col = c("red", "black"), pch = 1, bty = "n")
    legend("topright", legend=c("Mean", "Median"), col = "black", pch = c(1,2), bty = "n", inset = c(.15,.09))
}
png(file = paste0("/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/median_migration_curve_zoom_", name_append,"_",dat_rem, ".png"), res = 800, units = "in", height = 7, width = 7)
plotot()
dev.off()




# plot mean demography
pops_to_exclude <- "JEJ-C119|KoeJ-C212|VAJ-C24|HHJ-C162|SBJ-C52|LVJ-C76|ULJ-C122|SKJ-C16" #list them, separate by |
pops_to_exclude_later <- "ULJ-C122|SKJ-C16"
demo_files <- list.files(path = "/faststorage/project/EcoGenetics/people/Jesper_Bechsgaard/Entomobrya_nicoleti_paper1/figures/", pattern = "nosingletons.final.summary", full.names = TRUE)
demo_files_double <- demo_files[grep(paste0(pops_to_exclude), demo_files, invert = TRUE)]
demo_dat_double <- lapply(demo_files_double, read.table, header = TRUE)
demo_dat_names_double <- lapply(demo_files_double, function(x) sapply(strsplit(basename(x), "_"), "[[", 2))
if (any(grepl("ae", demo_dat_names_double, fixed = T))){       demo_dat_names_double <- gsub("ae", "A", demo_dat_names_double)    }
if (any(grepl("aa", demo_dat_names_double,  fixed = T))){       demo_dat_names_double <- gsub("aa", "A", demo_dat_names_double)    }
if (any(grepl("oe", demo_dat_names_double, fixed = T))){       demo_dat_names_double <- gsub("oe", "O", demo_dat_names_double)    }
demo_dat_double <- demo_dat_double[grep("SKJ",demo_dat_names_double )]
demo_dat_names_double <- demo_dat_names_double[grep("SKJ",demo_dat_names_double )]
names(demo_dat_double) <- demo_dat_names_double
str(demo_dat_double)
demo_dat_double <- lapply(demo_dat_double, function(x){ colnames(x)[which(colnames(x) == "year")] <- "K_year"; return(x)})
demo_dat_double <- lapply(demo_dat_double, function(x){ x[,which(colnames(x) == "K_year")] <- x[,which(colnames(x) == "K_year")]/1000; return(x)})

demo_files <- demo_files[grep(paste0(pops_to_exclude,"|",pops_to_exclude_later), demo_files, invert = TRUE)]
demo_dat <- lapply(demo_files, read.table, header = TRUE)

demo_dat_names <- lapply(demo_files, function(x) sapply(strsplit(basename(x), "-|_"), "[[", 2))
if (any(grepl("ae", demo_dat_names, fixed = T))){       demo_dat_names <- gsub("ae", "A", demo_dat_names)    }
if (any(grepl("aa", demo_dat_names,  fixed = T))){       demo_dat_names <- gsub("aa", "A", demo_dat_names)    }
if (any(grepl("oe", demo_dat_names, fixed = T))){       demo_dat_names <- gsub("oe", "O", demo_dat_names)    }
all( table(demo_dat_names) == 1 )
names(demo_dat) <- demo_dat_names
lapply(demo_dat, tail)
demo_dat <- lapply(demo_dat, function(x){ colnames(x)[which(colnames(x) == "year")] <- "K_year"; return(x)})
demo_dat <- lapply(demo_dat, function(x){ x[,which(colnames(x) == "K_year")] <- x[,which(colnames(x) == "K_year")]/1000; return(x)})
demo_dat <- lapply(demo_dat, function(x){ x[,which(colnames(x) == "Ne_12.5.")] <- x[,which(colnames(x) == "Ne_12.5.")]/100000; return(x)})
demo_dat <- lapply(demo_dat, function(x){ x[,which(colnames(x) == "Ne_87.5.")] <- x[,which(colnames(x) == "Ne_87.5.")]/100000; return(x)})
demo_dat <- lapply(demo_dat, function(x){ x[,which(colnames(x) == "Ne_median")] <- x[,which(colnames(x) == "Ne_median")]/100000; return(x)})

str(demo_dat)

for( element in seq(1, NROW(demo_dat))){
    if( element == 1){
        data_Ne <- demo_dat[[element]]$Ne_median
        ne_dat <- data.frame(demo_dat[[element]]$K_year, data_Ne)
        colnames(ne_dat)[1] <- "K_year"
        colnames(ne_dat)[element+1] <- paste0(names(demo_dat)[element],"_Ne")
    }else{
        data_Ne <- demo_dat[[element]]$Ne_median
        ne_dat <- cbind(ne_dat, data_Ne)
        colnames(ne_dat)[element+1] <- paste0(names(demo_dat)[element],"_Ne")
    }
}

z_trans_dat <- scale(ne_dat[,-1])
mean_traj <- rowMeans(ne_dat[,-1])

ylims = c(min(c(z_trans_dat, mean_traj)), max(c(z_trans_dat, mean_traj)))
#ylims = c(min(z_trans_dat), max(z_trans_dat))
library(scales)
# plot Ne trajectories with both mean of population Medians, and z-transformed trajectories (using scale) 
png(paste0("/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/Ne_stairway_trajectory_z_and_mean", name_append, ".png"), res = 600, height = 7, width = 7, units = "in")
plot(ne_dat[,1], z_trans_dat[,1], type = 'l', col = "grey", ylim = ylims, ylab = "z-transformed Ne (and blue, Nex 10^5)", xlab = "Thousand years bp")
for(colu in seq(2,ncol(z_trans_dat))){
    lines(ne_dat[,1], z_trans_dat[,colu], col = "grey", ylim = ylims)
}
lines(ne_dat[,1], rowMeans(z_trans_dat), ylim = ylims, lwd = 5)
lines(ne_dat[,1], mean_traj, ylim = ylims, lwd = 5, col = "blue")
dev.off()
 # plot only z transformed
ylims = c(min(z_trans_dat), max(z_trans_dat))
png(paste0("/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/Ne_stairway_trajectory_z_trans", name_append, ".png"), res = 600, height = 7, width = 7, units = "in")
plot(ne_dat[,1],z_trans_dat[,1], type = 'l', col = alpha("orange", .3), ylim = ylims, ylab = "z-transformed Ne", xlab = "Thousand years bp")
for(colu in seq(2,ncol(z_trans_dat))){
    lines(ne_dat[,1], z_trans_dat[,colu], col = alpha("orange", .3), ylim = ylims)
}
lines(ne_dat[,1], rowMeans(z_trans_dat), ylim = ylims, lwd = 5, col = "darkorange")
dev.off()



# landscape
library(viridis)
data <- as.data.frame(readRDS(file = "/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/Entomobrya_nicoleti/EntNic/figures/landscape/DF1.rds"))
head(data)
data$Km2 <- data$Km2/43094*100
data$Year <- abs(data$Year-2000)
head(data)
col_dat<- data.frame(unique(data$Landuse), rev(c(viridis(4))))

png("/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/landscape.png", res = 600, height = 7, width = 7, units = "in")
par(lwd = 1.7)
xlimis = c(0,10000)
plot(data[which(data$Landuse == col_dat[1,1]),"Year"], data[which(data$Landuse == col_dat[1,1]),"Km2"], type = "l", col = col_dat[1,2], xlab = NA, ylab = NA, xlim = xlimis)
lines(data[which(data$Landuse == col_dat[2,1]),"Year"], data[which(data$Landuse == col_dat[2,1]),"Km2"], col = col_dat[2,2], xlim = xlimis)
lines(data[which(data$Landuse == col_dat[3,1]),"Year"], data[which(data$Landuse == col_dat[3,1]),"Km2"], col = col_dat[3,2], xlim = xlimis)
lines(data[which(data$Landuse == col_dat[4,1]),"Year"], data[which(data$Landuse == col_dat[4,1]),"Km2"], col = col_dat[4,2], xlim = xlimis)
title(ylab = "Percentage of Landcover across Denmark", xlab = "Years bp")
#legend("topleft", legend = col_dat[,1], col = col_dat[,2], lty = 1)
legend("topright", legend = col_dat[,1], col = col_dat[,2], lty = 1, bty = "n")
dev.off()








# try to combine in one graph?
    # Agri, stairway (Ne), migration
barplot(rep(1,5), col = inferno(5, begin = .42, end = .62)[c(1,5)])
dev.off()
col2rgb(inferno(5, begin = .42, end = .62)[c(1,5)])
#[1] "#AD305DFF" "#D94D3DFF"
#       [,1] [,2]
# red    155  226
# green   41   87
# blue   100   52

col_pan <- rev(inferno(3, begin = .2, end = .8))
plotot <- function(){
    par(mfrow = c(1,1), mar = c(1.5,5,0,0), oma = c(6, 4, 1, 1))
    xlimis = c(0,4000);     xpd_is = FALSE;     mgp_is = c(3,1,0) # margin line for ax title, labels and line (3,1,0)
    cex_adj = .9
    # agriculture
    par(lwd = 2)
    crop_dat <- data[which(data$Landuse == "cropland"),]
    crop_dat <- crop_dat[order(crop_dat$Year),]
    plot(crop_dat$Year, crop_dat$Km2, type = "l", 
        col = "white", xlab = NA, ylab = NA, xlim = xlimis, xpd = xpd_is, xaxt = 'n', yaxt = 'n', frame.plot = F, 
        ylim = range(pretty(c(0, crop_dat$Km2))))
    #plot(data[which(data$Landuse == "cropland"),"Year"], data[which(data$Landuse == "cropland"),"Km2"], type = "l", 
    #   col = "white", xlab = NA, ylab = NA, xlim = xlimis, xpd = xpd_is, xaxt = 'n', yaxt = 'n', frame.plot = F, 
    #  ylim = range(pretty(c(0, data[which(data$Landuse == "cropland"),"Km2"]))))
    axis(2, las = 2, line = 0, col = col_pan[1], cex.axis = cex_adj)
    #rect( 7000/3, -10, 9500/3, 80, col = "grey90", border = NULL, lwd = 0, xpd = F)
    #rect( 100, -10, 650, 80, col = "grey90", border = NULL, lwd = 0, xpd = F)
    lines(crop_dat$Year, crop_dat$Km2, col = col_pan[1], xlim = xlimis, xpd = xpd_is)
    #lines(data[which(data$Landuse == "cropland"),"Year"], data[which(data$Landuse == "cropland"),"Km2"], col = col_pan[1], xlim = xlimis, xpd = xpd_is)

    # migration
    median_frame <- aggregate(df_sort, list(df_sort$migration_divide), FUN = median, rm.na = T)
    median_frame$mig_anc <- median_frame$mig_anc*100000
    median_frame$mig_rec <- median_frame$mig_rec*100000
    migmax <- max(c(median_frame$mig_rec, median_frame$mig_anc))
    par( new = T)
    plot(median_frame$Group.1, median_frame$mig_anc, #ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), migmax), 
        col = col_pan[2], ylab = NA, xlab = NA, frame.plot = F, 
        xlim = xlimis, xpd = xpd_is, type = "l", xaxt = 'n', yaxt = 'n', ylim = range(pretty(c(0, migmax))) )
    axis(2, las = 2, line = 3, col = col_pan[2], cex.axis = cex_adj)
    #mtext(expression("Gene flow" ~ (10^-5)), side = 2, adj = -.7, padj = -3.5, xpd = NA)    #title(ylab = expression("Migration rate" ~ (10^-5)), xpd = NA, mgp = mgp_is, line = 6.5)


    lines(median_frame$Group.1, median_frame$mig_rec, ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), migmax), 
        col = col_pan[2], xlim = xlimis, xpd = xpd_is)
    points(median_frame$Group.1, median_frame$mig_rec, ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), migmax), 
        col = col_pan[2], xlim = xlimis, xpd = xpd_is, cex = .4, pch = 19)

    # Ne trajectories (mean of medians)
    year_vec <- ne_dat[,1]*1000
    year_vec <- year_vec[which(year_vec < 4500)]
    ylims = range(pretty(c(min(rowMeans(z_trans_dat[1:NROW(year_vec),])), max(rowMeans(z_trans_dat[1:NROW(year_vec),])))))

    par( new = T)
    plot(year_vec, rowMeans(z_trans_dat[1:NROW(year_vec),]), col = col_pan[3], xpd = xpd_is, xlim = xlimis, type = 'l', yaxt = 'n',
        xlab = NA, ylab = NA, frame.plot = F, ylim = ylims)
    axis(2, las = 2, line = 6, col = col_pan[3], cex.axis = cex_adj)
    title(xlab = "Years bp", xpd = NA, mgp = c(2.5,1,0))
    #mtext("        Ne\n(z-transformed)", side = 2, adj = -.65, padj = -4, xpd = NA)    #title(ylab = "Ne (z-transformed)", xpd = NA, mgp = mgp_is, line = 9.5)

    text(-330, -.78, labels=c("Percent\nCropland cover"), col = col_pan[1], xpd = NA, srt = 90, adj = 0)
    text(-600, -.78, labels=c(expression("Gene flow" ~ (10^-5))), col = col_pan[2], xpd = NA, srt = 90, adj = 0)
    text(-900, -.78, labels=c(expression("N"[e] ~ "(z-transformed)")), col = col_pan[3], xpd = NA, srt = 90, adj = 0)

    legend(900, y=-.63, legend=c( "Ancestral gene flow", "Recent gene flow", "Cropland cover", expression("N"[e] ~ "trajectory")),
        col = col_pan[c(2,2,1,3)], lty = 1, pch = c(NA, 19, NA, NA), bty = "n", xpd = NA, ncol = 2,
        xjust = 0, text.width = 1000)
}
png(file = paste0("/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/median_migration_Ne_agri_ONE1_", name_append,"_",dat_rem,".png"), res = 800, units = "in", height = 6, width = 10)
plotot()
dev.off()

#inset = -.33





# simplify with fitted lines

col_pan <- rev(inferno(3, begin = .2, end = .8))
plotot <- function(){
    par(mfrow = c(1,1), mar = c(1.5,5,0,0), oma = c(6, 4, 1, 1))
    xlimis = c(0,4000);     xpd_is = FALSE;     mgp_is = c(3,1,0) # margin line for ax title, labels and line (3,1,0)
    cex_adj = .9
    # agriculture
    par(lwd = 2)
    crop_dat <- data[which(data$Landuse == "cropland"),]
    crop_dat <- crop_dat[order(crop_dat$Year),]
    plot(crop_dat$Year, crop_dat$Km2, type = "l", 
        col = "white", xlab = NA, ylab = NA, xlim = xlimis, xpd = xpd_is, xaxt = 'n', yaxt = 'n', frame.plot = F, 
        ylim = range(pretty(c(0, crop_dat$Km2))))
    #plot(data[which(data$Landuse == "cropland"),"Year"], data[which(data$Landuse == "cropland"),"Km2"], type = "l", 
    #   col = "white", xlab = NA, ylab = NA, xlim = xlimis, xpd = xpd_is, xaxt = 'n', yaxt = 'n', frame.plot = F, 
    #  ylim = range(pretty(c(0, data[which(data$Landuse == "cropland"),"Km2"]))))
    axis(2, las = 2, line = 0, col = col_pan[1], cex.axis = cex_adj)
    #rect( 7000/3, -10, 9500/3, 80, col = "grey90", border = NULL, lwd = 0, xpd = F)
    #text(7000/3 + ((9500/3-7000/3)/2), 70, labels = "Beginning of agriculture", cex = .7, col = "grey50")
    #rect( 50, -10, 650, 80, col = "grey90", border = NULL, lwd = 0, xpd = F)
    #text( 50 + ((650-50)/2), 70, labels = "Increasing agriculture", cex = .7, col = "grey50")
    # crop dat:
    lines(crop_dat$Year, crop_dat$Km2, col = col_pan[1], xlim = xlimis, xpd = xpd_is)
    #lines(data[which(data$Landuse == "cropland"),"Year"], data[which(data$Landuse == "cropland"),"Km2"], col = col_pan[1], xlim = xlimis, xpd = xpd_is)

    # migration
    median_frame <- aggregate(df_sort, list(df_sort$migration_divide), FUN = median, rm.na = T)
    median_frame$mig_anc <- median_frame$mig_anc*100000
    median_frame$mig_rec <- median_frame$mig_rec*100000
    migmax <- max(c(median_frame$mig_rec, median_frame$mig_anc))
    par( new = T)
    #plot(median_frame$Group.1, median_frame$mig_anc, #ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), migmax), 
    plot(NULL, NULL, #ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), migmax), 
        col = col_pan[2], ylab = NA, xlab = NA, frame.plot = F, 
        xlim = xlimis, xpd = xpd_is, type = "l", xaxt = 'n', yaxt = 'n', ylim = range(pretty(c(0, migmax))) )


    #}
    # recent:
    func_values <- predict(smooth.spline(x=median_frame$Group.1, y=median_frame$mig_rec), type = "response")
    func_values <- func_values$y
    print(func_values)
    z_vec <- vector()
    count_vec <- vector()
    sumo_vec <- vector()
    last_val_vec <- vector()
    for(each in seq(NROW(func_values))){
        if( each > 2){
            last_val <- func_values[each]
            # last_val = (sum(z_vec) + z_val)/(NROW(z_vec) + 1)
            # last_val = (sum(z_vec) + z_val)/(NROW(z_vec) + 1)
            # last_val * NROW(z_vec + 1) = sum(z_vec ) + z_val
            # last_val * NROW(z_vec + 1) - sum(z_vec ) =  z_val
            sum_o <- sum(z_vec[1:(each-1)])
            count_o <- NROW(z_vec[1:(each)])
            # avg_n <- sum(new_values[1:(each-1)]+last_val)/NROW(new_values[1:(each)])
            z_val <- (last_val * count_o) - sum_o
            z_vec <- append(z_vec, z_val)
            count_vec <- append(count_vec, count_o)
            sumo_vec <- append(sumo_vec, sum_o)
            last_val_vec <- append(last_val_vec, last_val)
        }else{
            z_vec <- append(z_vec, func_values[each])
        }
    }

    print(NROW(z_vec)); print(NROW(median_frame$Group.1))
    lines(median_frame$Group.1, (z_vec), col = col_pan[2], lty = 1, ylim = ylims, xlim = xlimis)
    lines((median_frame$Group.1-800)[-which((median_frame$Group.1-800) < 0)], (z_vec)[-which((median_frame$Group.1-800) < 0)], col = col_pan[2], lty = 2, ylim = ylims, xlim = xlimis, xpd = F)
    #lines(median_frame$Group.1, func_values, col = col_pan[2], lty = 3, ylim = ylims, xlim = xlimis)

    # calculate for polygon
    x <- c(median_frame$Group.1, rev(median_frame$Group.1-800))
    y <- c((z_vec), rev(z_vec))
    polygon(x,y, col = alpha(col_pan[2], .1), border = NA)
    arrows(x0 = 2600-10, x1 = 2200-10, y0 = .23, y1 = .23, length = .08,  col = alpha(col_pan[2], .7), lwd = .8)
    text( 2400-10, .26, labels = "X years", col = col_pan[2], cex = .7)
    #print((z_vec)/func_values)


    axis(2, las = 2, line = 3, col = col_pan[2], cex.axis = cex_adj)

    # Ne trajectories (mean of medians)
    year_vec <- ne_dat[,1]*1000
    year_vec <- year_vec[which(year_vec < 4500)]
    ylims = range(pretty(c(min(rowMeans(z_trans_dat[1:NROW(year_vec),])), max(rowMeans(z_trans_dat[1:NROW(year_vec),])))))

    par( new = T)
    plot(year_vec, rowMeans(z_trans_dat[1:NROW(year_vec),]), col = col_pan[3], xpd = xpd_is, xlim = xlimis, type = 'l', yaxt = 'n',
        xlab = NA, ylab = NA, frame.plot = F, ylim = ylims)
    axis(2, las = 2, line = 6, col = col_pan[3], cex.axis = cex_adj)
    title(xlab = "Years bp", xpd = NA, mgp = c(2.5,1,0))
    #mtext("        Ne\n(z-transformed)", side = 2, adj = -.65, padj = -4, xpd = NA)    #title(ylab = "Ne (z-transformed)", xpd = NA, mgp = mgp_is, line = 9.5)

    text(-255, -.8, labels=c("Cropland cover (%)"), col = col_pan[1], xpd = NA, srt = 90, adj = 0)
    text(-585, -.8, labels=c(expression("Gene flow" ~ (10^-5))), col = col_pan[2], xpd = NA, srt = 90, adj = 0)
    text(-900, -.8, labels=c(expression("N"[e] ~ "(z-transformed)")), col = col_pan[3], xpd = NA, srt = 90, adj = 0)

    #legend(900, y=-.63, legend=c( "Ancestral gene flow", "Recent gene flow", "Cropland cover", "Ne trajectory"),
    #    col = col_pan[c(2,2,1,3)], lty = 1, pch = c(NA, 19, NA, NA), bty = "n", xpd = NA, ncol = 2,
    #    xjust = 0, text.width = 1000)
    legend(900, y=-.62, legend=c( "Predicted gene flow", "Corrected gene flow", "Cropland cover", expression("N"[e] ~ "trajectory")),
        col = col_pan[c(2,2,1,3)], lty = c(1,2,1,1), bty = "n", xpd = NA, ncol = 2,
        xjust = 0, text.width = 1000)
}
png(file = paste0("/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/median_migration_Ne_agri_ONE_FITline_simple_polygon_", name_append,"_",dat_rem,".png"), res = 800, units = "in", height = 6, width = 10)
#png(file = paste0("/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/median_migration_Ne_agri_ONE_FITline_simple_poly_", name_append,"_",dat_rem,".png"), res = 800, units = "in", height = 6, width = 10)
#png(file = paste0("/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/median_migration_Ne_agri_ONE_FITline_simple_", name_append,"_",dat_rem,".png"), res = 800, units = "in", height = 6, width = 10)
plotot()
dev.off()

                                                              


#panel for supplement












###############################
### get the simulated SFS runs:
###############################

sim_list_of_files <- list.files(pattern="bestlikelyhoods" , path="/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/simulated_SFS/", recursive = TRUE, full.names = TRUE)
NROW(sim_list_of_files)

    # only files from one pop comparison
#file <- "/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/2dSFS/2dSFS_EntNic_aaRJ_C225_vs_EntNic_aeRoe-C36_1_vs_2//900/2dSFS_EntNic_aaRJ_C225_vs_EntNic_aeRoe-C36_1_vs_2_900migDiv_sfs2d/2dSFS_EntNic_aaRJ_C225_vs_EntNic_aeRoe-C36_1_vs_2_900migDiv_sfs2d.bestlhoods"
sim_df <- data.frame()
for(file in sim_list_of_files){
    datas <- read.table(file, header = T)
    name <- basename(file)
    count <- NROW(unlist(strsplit(name, "_")))
    migration_divide <- unlist(strsplit(name, "_"))[count-5]
    migration_divide <- as.integer(sapply(strsplit(migration_divide, "mig"), "[[", 1))
    #migration_divide <- migration_divide/3
    # get pop pair
    pop_pair <- gsub("2dSFS_", "", name)
    pop_pair <- gsub("_[0-9]*_vs_[0-9]*_[0-9]*migDiv_sfs2d.bestlikelyhoods", "", pop_pair)
    
    # make wide format
    datas1 <- data.frame(pop_pair, migration_divide, datas)
    sim_df <- rbind(sim_df, datas1)
}


str(sim_df)
sim_df_sort <- sim_df[order(sim_df[,1], sim_df[,2], method = "auto"),]
nrow(sim_df_sort)
# 1] 4166
sim_df_sort_alldat <- sim_df_sort
sim_df_sort <- sim_df_sort[which(sim_df_sort$migration_divide <= dat_rem),]

# remove pairs with split time below 12000
totl_nrow <- nrow(sim_df_sort)
print(totl_nrow)
sim_df_sort <- sim_df_sort[which(sim_df_sort$TDIV > dat_rem),]
nrow(sim_df_sort)
#[1] 3772
nrow(sim_df_sort)/totl_nrow
# [1] 0.9054249

# translate into years
sim_df_sort$migration_divide <- sim_df_sort$migration_divide/3

# add simulation type to df
sim_df_sort$sim_type <- sapply(strsplit(sim_df_sort$pop_pair, "\\.|rep_[0-9]_|rep_[0-9][0-9]_"), "[[", 2)
sim_df_sort_ne <- sim_df_sort[,c("migration_divide", "NPOP1", "NPOP2", "NANC", "sim_type")]
sim_df_sort <- sim_df_sort[,c("migration_divide", "mig_rec", "mig_anc", "sim_type")]
# real data
mean_frame <- aggregate(df_sort, list(df_sort$migration_divide), FUN = mean, rm.na = T)
median_frame <- aggregate(df_sort, list(df_sort$migration_divide), FUN = median, rm.na = T)
table (df_sort$migration_divide)
# constant_migration
constant_migration_df <- sim_df_sort[ which(sim_df_sort$sim_type == "constant_migration" ),] 
constant_migration_mean_frame <- aggregate(constant_migration_df, list(constant_migration_df$migration_divide), FUN = mean, rm.na = T)
constant_migration_median_frame <- aggregate(constant_migration_df, list(constant_migration_df$migration_divide), FUN = median, rm.na = T)
table(constant_migration_df$migration_divide)
# no_migration
no_migration_df <- sim_df_sort[ which(sim_df_sort$sim_type == "no_migration" ),] 
no_migration_mean_frame <- aggregate(no_migration_df, list(no_migration_df$migration_divide), FUN = mean, rm.na = T)
no_migration_median_frame <- aggregate(no_migration_df, list(no_migration_df$migration_divide), FUN = median, rm.na = T)
table(no_migration_df$migration_divide)
# decreasing_migration
decreasing_migration_df <- sim_df_sort[ which(sim_df_sort$sim_type == "decreasing_migration" ),] 
decreasing_migration_mean_frame <- aggregate(decreasing_migration_df, list(decreasing_migration_df$migration_divide), FUN = mean, rm.na = T)
decreasing_migration_median_frame <- aggregate(decreasing_migration_df, list(decreasing_migration_df$migration_divide), FUN = median, rm.na = T)
table(decreasing_migration_df$migration_divide)

ylims_mean = c( min(c(mean_frame$mig_rec, mean_frame$mig_anc, constant_migration_mean_frame$mig_rec, constant_migration_mean_frame$mig_anc, no_migration_mean_frame$mig_rec, no_migration_mean_frame$mig_anc, decreasing_migration_mean_frame$mig_rec, decreasing_migration_mean_frame$mig_anc)),
            max(c(mean_frame$mig_rec, mean_frame$mig_anc, constant_migration_mean_frame$mig_rec, constant_migration_mean_frame$mig_anc, no_migration_mean_frame$mig_rec, no_migration_mean_frame$mig_anc, decreasing_migration_mean_frame$mig_rec, decreasing_migration_mean_frame$mig_anc)))
ylims_median = c( min(c(median_frame$mig_rec, median_frame$mig_anc, constant_migration_median_frame$mig_rec, constant_migration_median_frame$mig_anc, no_migration_median_frame$mig_rec, no_migration_median_frame$mig_anc, decreasing_migration_median_frame$mig_rec, decreasing_migration_median_frame$mig_anc)),
            max(c(median_frame$mig_rec, median_frame$mig_anc, constant_migration_median_frame$mig_rec, constant_migration_median_frame$mig_anc, no_migration_median_frame$mig_rec, no_migration_median_frame$mig_anc, decreasing_migration_median_frame$mig_rec, decreasing_migration_median_frame$mig_anc)))
ylims <- c(min(ylims_mean, ylims_median), max(ylims_mean, ylims_median))
ylims <- c(ylims[1], .00004)
ylims <- c(ylims[1], .00006)
colors_are <- c("red", "black")
colors_are <- plasma(4, end=.8)
pch_are <- c(4, 1)



xlima = c(0,14000/3)
# in three graphs
plotthree <- function(){
    ylims = c(0,6*10^-5)
    colors_are <- plasma(4, end=.8)
    par(mfrow = c(3,1), mar = c(1,5.1,1,1), oma = c(4,0,3,0))
    dat_list <- list(median_frame, constant_migration_median_frame, decreasing_migration_median_frame, no_migration_median_frame)
    names(dat_list) <- c("median_frame", "constant_migration_median_frame", "decreasing_migration_median_frame", "no_migration_median_frame")
    names(dat_list) <- paste0("simulated_", names(dat_list))
    names(dat_list)[1] <- "observed_data"
    names(dat_list) <- gsub("migration", "gene flow", gsub("_", " ", gsub("_median_frame", "", names(dat_list))))
    substr(names(dat_list), 1, 1) <- toupper(substr(names(dat_list), 1, 1)) # change first letter to capital
    # Constant
    plot(NULL, xlim = xlima, ylim = ylims, ylab = "Gene flow rate", xlab = NULL, xaxt = 'n', las = 2, mgp = c(4,1,0))
    sequence_is <- seq(1, NROW(dat_list))
    #for( number_is in sequence_is[c(1,2)] ){
    for( number_is in sequence_is[c(2)] ){
        migration_frame <- dat_list[[number_is]]
        points(migration_frame$Group.1, migration_frame$mig_anc, col = alpha(colors_are[number_is],.7), pch = pch_are[1], xlim = xlima, ylim = ylims)
        points(migration_frame$Group.1, migration_frame$mig_rec, col = alpha(colors_are[number_is],.7), pch = pch_are[2], xlim = xlima, ylim = ylims)
    }
    
    abline(h = 0.00006, col = "black", lty = 2)
    #abline(h = 0.00006, col = colors_are[2], lty = 4)
    title(main = names(dat_list)[2], col.main = colors_are[2], adj = 0.02, line = -2.5, xpd = NA)
    axis(1, labels = NA)
    mtext("a", side = 3, adj = -.11, line = .5, cex = 1.1)
    #axis(1, labels = NA, col = "grey")
    #box(col = "grey")
    #legend("top", legend=c("Ancestral gene flow", "Recent gene flow", names(dat_list)[c(1)], "Gene flow rate in simulation"), 
    legend("top", legend=c(expression("Gene flow t"[1] ~ "(historical)"), expression("Gene flow t"[2] ~ "(recent)"), "Gene flow rate in simulation"), 
        #col = alpha(c(rep(colors_are[2],2), colors_are[1], colors_are[2]), .7),   
        col = "black",   
        pch = c(pch_are, NA),
        lty = c(rep(NA, 2), 2), 
        text.col = "black", # text.col = c(rep("black",2), colors_are[1], "black"),
        bty = "n", xpd = NA, inset = c(0, -.2), ncol = 2)
    #legend("topleft", legend=c(names(dat_list)[c(1,2)], "Gene flow rate in simulation"), col = colors_are[c(1,2,2)], 
    #    pch = c(pch_are[2], pch_are[2], NA), lty = c(NA,NA, 1), bty = "n", xpd = NA, inset = c(0, .03))

    # Decreasing
    ylims = c(0,5*10^-5)
    plot(NULL, xlim = xlima, ylim = ylims, ylab = "Gene flow rate", xlab = NULL, xaxt = 'n', las = 2, mgp = c(4,1,0))
    #for( number_is in sequence_is[c(1,3)] ){
    for( number_is in sequence_is[c(3)] ){
        migration_frame <- dat_list[[number_is]]
        points(migration_frame$Group.1, migration_frame$mig_anc, col = alpha(colors_are[number_is],.7), pch = pch_are[1], xlim = xlima, ylim = ylims)
        points(migration_frame$Group.1, migration_frame$mig_rec, col = alpha(colors_are[number_is],.7), pch = pch_are[2], xlim = xlima, ylim = ylims)
    }
    segments(   x0 = c(0, 600, 3000, 6000, 9000, 12000, 15000)/3, 
                x1 = c(600, 3000, 6000, 9000, 12000, 15000, 45000)/3,
                y0 = c(0, 0.00001, 0.00002, 0.00003, 0.00004, 0.00005, 0.00006 ),
                y1 = c(0, 0.00001, 0.00002, 0.00003, 0.00004, 0.00005, 0.00006 ), 
                col = "black", lty = 2) 
                #col = colors_are[3], lty = 4) 
    axis(1, labels = NA)
    mtext("b", side = 3, adj = -.11, line = 0, cex = 1.1)
    #func_values <- predict(smooth.spline(x=migration_frame$Group.1, y=migration_frame$mig_rec), type = "response")
    #func_values <- func_values$y
    ##print(func_values)
    #z_vec <- vector()
    #count_vec <- vector()
    #sumo_vec <- vector()
    #last_val_vec <- vector()
    #for(each in seq(NROW(func_values))){
     #   if( each > 2){
      #      last_val <- func_values[each]
       #     # last_val = (sum(z_vec) + z_val)/(NROW(z_vec) + 1)
        #    # last_val = (sum(z_vec) + z_val)/(NROW(z_vec) + 1)
         #   # last_val * NROW(z_vec + 1) = sum(z_vec ) + z_val
          #  # last_val * NROW(z_vec + 1) - sum(z_vec ) =  z_val
           # sum_o <- sum(z_vec[1:(each-1)])
            #count_o <- NROW(z_vec[1:(each)])
#            # avg_n <- sum(new_values[1:(each-1)]+last_val)/NROW(new_values[1:(each)])
 #           z_val <- (last_val * count_o) - sum_o
  #          z_vec <- append(z_vec, z_val)
   #         count_vec <- append(count_vec, count_o)
    #        sumo_vec <- append(sumo_vec, sum_o)
     #       last_val_vec <- append(last_val_vec, last_val)
      #  }else{
       #     z_vec <- append(z_vec, func_values[each])
        #}
    #}

    #print(NROW(z_vec)); print(NROW(median_frame$Group.1))
    #lines(median_frame$Group.1, (z_vec), col = colors_are[3], lty = 2, ylim = ylims, xlim = xlima)
    #lines(median_frame$Group.1, func_values, col = colors_are[3], lty = 1, ylim = ylims, xlim = xlima)
    #legend("topleft", inset = c(0, .15), legend = c("Predicted gene flow (t2)", "Fitted line (t2)"), lty = c(2,1), col = colors_are[3], bty = "n")
    
    title(main = names(dat_list)[3], col.main = colors_are[3], adj = .02, line = -2, xpd = NA)
    #legend("topleft", legend=c(names(dat_list)[c(1,3)], "Gene flow rate in simulation"), col = colors_are[c(1,3,3)], 
     #   pch = c(pch_are[2], pch_are[2], NA), lty = c(NA,NA, 1), bty = "n", xpd = NA)

    # None
    plot(NULL, xlim = xlima, ylim = ylims, ylab = "Gene flow rate", las = 2, xaxt = 'n', mgp = c(4,1,0))
    abline(h = 0, col = "black", lty = 2)
    #for( number_is in sequence_is[c(1,4)] ){
    for( number_is in sequence_is[c(4)] ){
        migration_frame <- dat_list[[number_is]]
        points(migration_frame$Group.1, migration_frame$mig_anc, col = alpha(colors_are[number_is],.65), pch = pch_are[1], xlim = xlima, ylim = ylims)
        points(migration_frame$Group.1, migration_frame$mig_rec, col = alpha(colors_are[number_is],.65), pch = pch_are[2], xlim = xlima, ylim = ylims)
    }
    axis(1, xpd = NA)   
    mtext("c", side = 3, adj = -.11, line = 0, cex = 1.1)
    #axis(1, xpd = NA, col = "grey")
    #box(col = "grey")
    title(main = names(dat_list)[4], col.main = colors_are[4], adj = .02, line = -2, xpd = NA)
    title(xlab = "Years bp", xpd = NA)
    #legend("topleft", legend=c(names(dat_list)[c(1,4)]), col = colors_are[c(1,4)], 
     #   pch = pch_are[2], bty = "n", xpd = NA)
}
png(file = paste0("/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/median_migration_curve_EntNic_SIM_3graph","_",dat_rem,".png"), 
        res = 800, units = "in", height = 8, width = 6)
plotthree()
dev.off()

### OBS change!
# remove observed data
# change simulated lines to black? some other color?
# make legend simpler



mean_it <- aggregate(sim_df_sort_ne, list(sim_df_sort_ne$migration_divide), FUN = mean, rm.na = T)
median_it <- aggregate(sim_df_sort_ne, list(sim_df_sort_ne$migration_divide), FUN = median, rm.na = T)

sim_df_sort_ne_con <- sim_df_sort_ne[ which(sim_df_sort_ne$sim_type == "constant_migration" ),] 
mean_it_con <- aggregate(sim_df_sort_ne_con, list(sim_df_sort_ne_con$migration_divide), FUN = mean, rm.na = T)
median_it_con <- aggregate(sim_df_sort_ne_con, list(sim_df_sort_ne_con$migration_divide), FUN = median, rm.na = T)

sim_df_sort_ne_no <- sim_df_sort_ne[ which(sim_df_sort_ne$sim_type == "no_migration" ),] 
mean_it_no <- aggregate(sim_df_sort_ne_no, list(sim_df_sort_ne_no$migration_divide), FUN = mean, rm.na = T)
median_it_no <- aggregate(sim_df_sort_ne_no, list(sim_df_sort_ne_no$migration_divide), FUN = median, rm.na = T)

sim_df_sort_ne_dec <- sim_df_sort_ne[ which(sim_df_sort_ne$sim_type == "decreasing_migration" ),] 
mean_it_dec <- aggregate(sim_df_sort_ne_dec, list(sim_df_sort_ne_dec$migration_divide), FUN = mean, rm.na = T)
median_it_dec <- aggregate(sim_df_sort_ne_dec, list(sim_df_sort_ne_dec$migration_divide), FUN = median, rm.na = T)


# plot recent actross all simulation types

plottit <- function(){
    ylima = c(min(c(mean_it$NPOP1, mean_it$NPOP2, 22000)), max(c(mean_it$NPOP1, mean_it$NPOP2)))
    par(mar = c(5.1,5.1,1,1))
    cols_are <- inferno(3, begin = .2, end = .7)
    alpha_is = .6
    plot(NULL, xlim = c(0,5000), ylim = ylima, ylab = expression("N"[e] ~ "(recent, haploid)"), las = 2, mgp = c(4,1,0), type = 'p', xlab = NA, xaxt = 'n')
    #points(mean_it$Group.1, mean_it$NPOP1, col = alpha(cols_are[1],alpha_is))
    #points(mean_it$Group.1, mean_it$NPOP2, col = alpha(cols_are[1],alpha_is))
    lines(mean_it$Group.1, mean_it$NPOP1, col = alpha(cols_are[1],alpha_is))
    lines(mean_it$Group.1, mean_it$NPOP2, col = alpha(cols_are[1],alpha_is))
    lines(median_it$Group.1, median_it$NPOP1, col = alpha(cols_are[2],1))
    lines(median_it$Group.1, median_it$NPOP2, col = alpha(cols_are[2],1))
    abline(h = 22500, col = cols_are[3])
    axis(1, las = 1)
    title(xlab = "Generations bp")
    legend("topright", legend = c("Mean", "Median", expression("N"[e] ~ "in simulation")), col = cols_are, lty = 1, bty = 'n')
}
png(file = paste0("/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/NE_rec_curve_EntNic_SIM","_",dat_rem,".png"), 
    res = 800, units = "in", height = 5.5, width = 8)
plottit()
dev.off()

# plot recent actross all types
plottit <- function(){
    cols <- (plasma(4, begin = .1, end = .85))
    ylima = c(min(c(mean_it$NPOP1, mean_it$NPOP2, 22000)), 
            max(c(mean_it$NPOP1, mean_it$NPOP2)))
    alpha_is = .6
    par(mar = c(5.1,5.1,1,1))
    plot(NULL, xlim = c(0,5000), ylim = ylima, ylab = expression("N"[e] ~ "(recent, haploid)"), las = 2, mgp = c(4,1,0), type = 'p', xlab = NA, xaxt = 'n')
    points(mean_it_no$Group.1, mean_it_no$NPOP1, col = alpha(cols[1],alpha_is))
    points(mean_it_no$Group.1, mean_it_no$NPOP2, col = alpha(cols[1],alpha_is))
    lines(smooth.spline(mean_it_no$Group.1, mean_it_no$NPOP2), col = alpha(cols[1],1))
    #lines(mean_it_no$Group.1, mean_it_no$NPOP1, col = alpha(cols[1],alpha_is))
    #lines(mean_it_no$Group.1, mean_it_no$NPOP2, col = alpha(cols[1],alpha_is))
    #lines(median_it_no$Group.1, median_it_no$NPOP1, col = alpha(cols[1],alpha_is))
    #lines(median_it_no$Group.1, median_it_no$NPOP2, col = alpha(cols[1],alpha_is))

    points(mean_it_con$Group.1, mean_it_con$NPOP1, col = alpha(cols[2],alpha_is))
    points(mean_it_con$Group.1, mean_it_con$NPOP2, col = alpha(cols[2],alpha_is))
    lines(smooth.spline(mean_it_con$Group.1, mean_it_con$NPOP2), col = alpha(cols[2],1))
    #lines(mean_it_con$Group.1, mean_it_con$NPOP1, col = alpha(cols[2],alpha_is))
    #lines(mean_it_con$Group.1, mean_it_con$NPOP2, col = alpha(cols[2],alpha_is))
    #lines(median_it_con$Group.1, median_it_con$NPOP1, col = alpha(cols[2],alpha_is))
    #lines(median_it_con$Group.1, median_it_con$NPOP2, col = alpha(cols[2],alpha_is))

    points(mean_it_dec$Group.1, mean_it_dec$NPOP1, col = alpha(cols[3],alpha_is))
    points(mean_it_dec$Group.1, mean_it_dec$NPOP2, col = alpha(cols[3],alpha_is))
    lines(smooth.spline(mean_it_dec$Group.1, mean_it_dec$NPOP2), col = alpha(cols[3],1))
    #lines(mean_it_dec$Group.1, mean_it_dec$NPOP1, col = alpha(cols[3],alpha_is))
    #lines(mean_it_dec$Group.1, mean_it_dec$NPOP2, col = alpha(cols[3],alpha_is))
    #lines(median_it_dec$Group.1, median_it_dec$NPOP1, col = alpha(cols[3],alpha_is))
    #lines(median_it_dec$Group.1, median_it_dec$NPOP2, col = alpha(cols[3],alpha_is))
    abline(h = 22500, col = cols[4])
    title(xlab = "Generations bp")
    legend("topright", legend = c("No gene flow", "Constant gene flow", "Decreasing gene flow", expression("N"[e] ~ "in simulation")), 
            col = c(cols), lty = 1, bty = "n")
    axis(1, las = 1)
}
png(file = paste0("/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/NE_rec_types_curve_EntNic_SIM","_",dat_rem,".png"), 
    res = 800, units = "in", height = 7, width = 8)
plottit()
dev.off()


# plot ancestral actross all simulation types
plottit <- function(){
    ylima = c(min(c(mean_it$NANC, median_it$NANC)), max(c(mean_it$NANC, median_it$NANC)))
    par(mar = c(5.1,5.1,1,1))
    alpha_is = .6
    plot(NULL, xlim = c(0,5000), ylim = ylima, ylab = expression("N"[e] ~ "(ancestral, haploid)"), las = 2, mgp = c(4,1,0), type = 'p', xlab = NA)
    points(mean_it$Group.1, mean_it$NANC, col = alpha("green4",alpha_is))
    lines(mean_it$Group.1, mean_it$NANC, col = alpha("green4",alpha_is))
    lines(median_it$Group.1, median_it$NANC, col = alpha("green4",.9))
    lines(median_it$Group.1, median_it$NANC, col = alpha("green4",.9))
    abline(h = 405000, col = "orange")
    title(xlab = "Generations bp")
}
png(file = paste0("/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/NE_anc_curve_EntNic_SIM","_",dat_rem,".png"), 
    res = 800, units = "in", height = 5, width = 7)
plottit()
dev.off()

# plot ancestral actross NO MIGRATION
plottit <- function(){
    cols <- plasma(3, begin = .2, end = .8)
    alpha_is = .6
    ylima = c(min(c(mean_it_no$NANC, mean_it_con$NANC, mean_it_dec$NANC)), max(c(mean_it_no$NANC, mean_it_con$NANC, mean_it_dec$NANC)))
    par(mar = c(5.1,5.1,1,1))
    plot(NULL, xlim = c(0,5000), ylim = ylima, ylab = expression("N"[e] ~ "(ancestral, haploid)"), las = 2, mgp = c(4,1,0), type = 'p', xlab = NA)
    points(mean_it_no$Group.1, mean_it_no$NANC, col = alpha(cols[1],alpha_is))
    lines(mean_it_no$Group.1, mean_it_no$NANC, col = alpha(cols[1],alpha_is))
    lines(median_it_no$Group.1, median_it_no$NANC, col = alpha(cols[1],alpha_is))

    points(mean_it_con$Group.1, mean_it_con$NANC, col = alpha(cols[2],alpha_is))
    lines(mean_it_con$Group.1, mean_it_con$NANC, col = alpha(cols[2],alpha_is))
    lines(median_it_con$Group.1, median_it_con$NANC, col = alpha(cols[2],alpha_is))

    points(mean_it_dec$Group.1, mean_it_dec$NANC, col = alpha(cols[3],alpha_is))
    lines(mean_it_dec$Group.1, mean_it_dec$NANC, col = alpha(cols[3],alpha_is))
    lines(median_it_dec$Group.1, median_it_dec$NANC, col = alpha(cols[3],alpha_is))

    abline(h = 405000, col = "blue")
    title(xlab = "Generations bp")
    legend("topright", legend = c("No gene flow", "Constant gene flow", "Decreasing gene flow", expression("Simulated recent N"[e])), col = c(cols, "blue"), lty = 1, bty = "n")
}
png(file = paste0("/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/NE_anc_types_curve_EntNic_SIM","_",dat_rem,".png"), 
    res = 800, units = "in", height = 7, width = 8)
plottit()
dev.off()












###############################
### get the simulated SFS runs: 3 rates!
###############################

sim_list_of_files_3rate <- list.files(pattern="bestlikelyhoods" , path="/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/simulated_SFS_3rates/", recursive = TRUE, full.names = TRUE)
NROW(sim_list_of_files_3rate)

    # only files from one pop comparison
#file <- "/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/2dSFS/2dSFS_EntNic_aaRJ_C225_vs_EntNic_aeRoe-C36_1_vs_2//900/2dSFS_EntNic_aaRJ_C225_vs_EntNic_aeRoe-C36_1_vs_2_900migDiv_sfs2d/2dSFS_EntNic_aaRJ_C225_vs_EntNic_aeRoe-C36_1_vs_2_900migDiv_sfs2d.bestlhoods"
sim_df_3rate <- data.frame()
for(file in sim_list_of_files_3rate){
    datas <- read.table(file, header = T)
    name <- basename(file)
    count <- NROW(unlist(strsplit(name, "_")))
    migration_divide <- unlist(strsplit(name, "_"))[count-5]
    migration_divide <- as.integer(sapply(strsplit(migration_divide, "mig"), "[[", 1))
    #migration_divide <- migration_divide/3
    # get pop pair
    pop_pair <- gsub("2dSFS_", "", name)
    pop_pair <- gsub("_[0-9]*_vs_[0-9]*_[0-9]*migDiv_sfs2d.bestlikelyhoods", "", pop_pair)
    
    # make wide format
    datas1 <- data.frame(pop_pair, migration_divide, datas)
    sim_df_3rate <- rbind(sim_df_3rate, datas1)
}


str(sim_df_3rate)
sim_df_3rate_sort <- sim_df_3rate[order(sim_df_3rate[,1], sim_df_3rate[,2], method = "auto"),]
nrow(sim_df_3rate_sort)
# 1] 4079
sim_df_sort_3rate_alldat <- sim_df_3rate_sort
sim_df_3rate_sort <- sim_df_3rate_sort[which(sim_df_3rate_sort$migration_divide <= dat_rem),]

# remove pairs with split time below 12000
totl_nrow <- nrow(sim_df_3rate_sort)
print(totl_nrow)
#[1] 1358
sim_df_3rate_sort <- sim_df_3rate_sort[which(sim_df_3rate_sort$TDIV > dat_rem),]
nrow(sim_df_3rate_sort)
#[1] 1316
nrow(sim_df_3rate_sort)/totl_nrow
# 0.9690722

# translate into years
sim_df_3rate_sort$migration_divide <- sim_df_3rate_sort$migration_divide/3

# add simulation type to df
sim_df_3rate_sort$sim_type <- sapply(strsplit(sim_df_3rate_sort$pop_pair, "\\.|rep_[0-9]_|rep_[0-9][0-9]_"), "[[", 2)
sim_df_3rate_sort_ne <- sim_df_3rate_sort[,c("migration_divide", "NPOP1", "NPOP2", "NANC", "sim_type")]
sim_df_3rate_sort <- sim_df_3rate_sort[,c("migration_divide", "mig_rec", "mig_int", "mig_anc", "sim_type")]
# real data
#mean_frame <- aggregate(sim_df_sort, list(sim_df_sort$migration_divide), FUN = mean, rm.na = T)
#median_frame <- aggregate(sim_df_sort, list(sim_df_sort$migration_divide), FUN = median, rm.na = T)
#table (sim_df_sort$migration_divide)
# decreasing_migration
decreasing_migration_df <- sim_df_3rate_sort[ which(sim_df_3rate_sort$sim_type == "decreasing_migration" ),] 
decreasing_migration_mean_frame <- aggregate(decreasing_migration_df, list(decreasing_migration_df$migration_divide), FUN = mean, rm.na = T)
decreasing_migration_median_frame <- aggregate(decreasing_migration_df, list(decreasing_migration_df$migration_divide), FUN = median, rm.na = T)
table(decreasing_migration_df$migration_divide)

ylims_mean = c( min(c(decreasing_migration_mean_frame$mig_rec, decreasing_migration_mean_frame$mig_anc, decreasing_migration_mean_frame$mig_int)),
            max(c(decreasing_migration_mean_frame$mig_rec, decreasing_migration_mean_frame$mig_anc, decreasing_migration_mean_frame$mig_int)))
ylims_median = c( min(c(decreasing_migration_median_frame$mig_rec, decreasing_migration_median_frame$mig_anc, decreasing_migration_median_frame$mig_int)),
            max(c(decreasing_migration_median_frame$mig_rec, decreasing_migration_median_frame$mig_anc, decreasing_migration_median_frame$mig_int)))
#ylims <- c(min(ylims_mean[,1], ylims_median[,1]), max(ylims_mean[,2], ylims_median[,2]))
#ylims <- c(min( ylims_median[,1]), max( ylims_median[,2]))
ylims <- c(ylims[1], .0002)
#ylims <- c(ylims[1], .00006)
colors_are <- c("red", "black")
colors_are <- plasma(4, end=.8)
pch_are <- c(4, 1)
plottit <- function(){
    #ylims = ylims_median
    ylims = ylims
    plot(NULL, xlim = c(0, 12000/3), ylim = ylims, ylab = "Gene flow rate", xlab = "Migration divide (Generations)")
    #for( number_is in seq(1, NROW(list(mean_frame, constant_migration_mean_frame, no_migration_mean_frame, decreasing_migration_mean_frame))) ){
    #   migration_frame <- list(mean_frame, constant_migration_mean_frame, no_migration_mean_frame, decreasing_migration_mean_frame)[[number_is]]
    #  print(number_is)
    # points(migration_frame$Group.1, migration_frame$mig_anc, col = alpha(colors_are[1],.3), pch = pch_are[number_is])
        #lines(migration_frame$Group.1, migration_frame$mig_anc, col = alpha(colors_are[1],.3))
    #    points(migration_frame$Group.1, migration_frame$mig_rec, col = alpha(colors_are[2], .3), pch = pch_are[number_is])
    #   lines(migration_frame$Group.1, migration_frame$mig_rec, col = alpha(colors_are[2], .3))
    #}
    dat_list <- list(decreasing_migration_median_frame)
    for( number_is in seq(1, NROW(dat_list)) ){
        migration_frame <- dat_list[[number_is]]
        print(number_is)
        points(migration_frame$Group.1, migration_frame$mig_anc, col = colors_are[number_is], pch = pch_are[1])
        #lines(migration_frame$Group.1, migration_frame$mig_anc, col = colors_are[number_is])
        lines(smooth.spline(migration_frame$Group.1, migration_frame$mig_anc), col = colors_are[number_is], lty = 3)
        points(migration_frame$Group.1, migration_frame$mig_rec, col = colors_are[number_is], pch = pch_are[2])
        #lines(migration_frame$Group.1, migration_frame$mig_rec, col = colors_are[number_is])
        lines(smooth.spline(migration_frame$Group.1, migration_frame$mig_rec), col = colors_are[number_is])
        points(migration_frame$Group.1, migration_frame$mig_int, col = colors_are[number_is], pch = 3)
        #lines(migration_frame$Group.1, migration_frame$mig_rec, col = colors_are[number_is])
        lines(smooth.spline(migration_frame$Group.1, migration_frame$mig_int), col = colors_are[number_is])
    }
    #abline(h = 0.00006, col = colors_are[2])
    segments(   x0 = c(0, 600, 3000, 6000, 9000, 12000, 15000)/3, 
                x1 = c(600, 3000, 6000, 9000, 12000, 15000, 45000)/3,
                y0 = c(0, 0.00001, 0.00002, 0.00003, 0.00004, 0.00005, 0.00006 ),
                y1 = c(0, 0.00001, 0.00002, 0.00003, 0.00004, 0.00005, 0.00006 ), 
                col = colors_are[3]) 
    legend("topleft", legend=c("Ancestral gene flow", "Recent gene flow", "Intermediate gene flow", "Simulated gene flow rate"), col = c(rep(colors_are[1],3), colors_are[3]), pch = c(pch_are, 3, NA), bty = "n", xpd = NA, lty = c(rep(NA, 3), 1))
    #legend("topleft", legend=c("Constant gene flow", "Decreasing gene flow"), lty = 1, col = colors_are[c(2,3)])
}
png(file = paste0( "/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/median_migration_curve_EntNic_SIM_3rates_",dat_rem,".png"), res = 800, units = "in", height = 7, width = 10)
plottit()
dev.off()





















###########################
### Panel on fastsimcoal for main
###########################
# model for estimaters
# gene flow dat




###########################
### Panel on fastsimcoal for supplement
###########################
# Ne data   -   Ne sim
# Split data    -   Split sim
# gene flow dat     -   model explaining x
    # gene flow without agri and ne, but with mig rec fitted


#df_sort_alldat
##### Ne data dat
library(tidyr)
head(df_sort_alldat)
cols_fsc <- (alpha(inferno(3, begin = .1, end = .38), .85))
df_sort_alldat$NPOP <- ((df_sort_alldat$NPOP1 + df_sort_alldat$NPOP2) / 2)
ne_data_plot <- function(){
    print(paste("Rec median, mean: ", median(df_sort_alldat[which(df_sort_alldat$TDIV > 14000),]$NPOP)/2, ", ", mean(df_sort_alldat[which(df_sort_alldat$TDIV > 14000),]$NPOP)/2))
    #[1] 43367.5
    #[1] 40543.04
    print(paste("Anc median, mean: ", median(df_sort_alldat[which(df_sort_alldat$TDIV > 14000),]$NANC)/2, ", ", mean(df_sort_alldat[which(df_sort_alldat$TDIV > 14000),]$NANC)/2))    
    print(paste("split median, mean: ", median(df_sort_alldat[which(df_sort_alldat$TDIV > 14000),]$TDIV)/2, ", ", mean(df_sort_alldat[which(df_sort_alldat$TDIV > 14000),]$TDIV)/2))    
    #[1] 2335525
    #[1] 2937668
    df_sort_alldat_NE_long <- gather(df_sort_alldat[which(df_sort_alldat$TDIV > 14000),], population, Ne, c(NPOP,NANC), factor_key=TRUE)
    print(nrow(df_sort_alldat[which(df_sort_alldat$TDIV > 14000),]))
    df_sort_alldat_NE_long$Ne <- (df_sort_alldat_NE_long$Ne/2)/100000
    # divide by 100000 to get on resonable scale
    # divide by 2 to get diploid (fastsimcoal estinates haploid pop sizes)
    boxplot(formula = Ne ~ population, data = df_sort_alldat_NE_long, col = cols_fsc[1:2], yaxt = "n", 
        names = c(expression(Avg(N[1] + N[2])), expression(N[Ancestral])), border = cols_fsc[1:2], medcol = "gray20", 
        ylab = expression("N"[e] ~ (10^5)), xlab = NULL)
    axis(2, las = 2)
}
#ne_dat()
ne_data_plot<- function(){
    print(paste("Rec median, mean: ", median(df_sort_alldat$NPOP/2), ", ", mean(df_sort_alldat$NPOP/2)))
    #[1] 43367.5
    #[1] 40543.04
    print(paste("Anc median, mean: ", median(df_sort_alldat$NANC/2), ", ", mean(df_sort_alldat$NANC/2)))  
    #[1] 2335525
    #[1] 2937668
    print(nrow(df_sort_alldat))
    df_sort_alldat_NE_long <- gather(df_sort_alldat, population, Ne, c(NPOP,NANC), factor_key=TRUE)
    df_sort_alldat_NE_long$Ne <- (df_sort_alldat_NE_long$Ne/2)/100000
    # divide by 100000 to get on resonable scale
    # divide by 2 to get diploid (fastsimcoal estinates haploid pop sizes)
    boxplot(formula = Ne ~ population, data = df_sort_alldat_NE_long, col = cols_fsc[1:2], yaxt = "n", 
        names = c(expression(Avg(N[1] + N[2])), expression(N[Ancestral])), border = cols_fsc[1:2], medcol = "gray20", 
        ylab = expression("N"[e] ~ (10^5)), xlab = NULL)
    axis(2, las = 2)
    
}
png(file = paste0("/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/panel_sup_fastsimcoal_part_a_ne_data.png"), 
    res = 800, units = "in", height = 8, width = 7)
par(mar= c(3.2, 4.6, 1, 1))
ne_data_plot()
dev.off()


##### Ne data SIM
library(tidyr)
head(sim_df_sort_alldat)
sim_df_sort_alldat$NPOP <- ((sim_df_sort_alldat$NPOP1 + sim_df_sort_alldat$NPOP2) / 2)
seg_col <- inferno(1, begin = .7, end = .72)
ne_sim <- function(){
    print(paste("Rec median, mean: ", median(sim_df_sort_alldat$NPOP), ", ", mean(sim_df_sort_alldat$NPOP)))    
    #[1] 41943.5
    #[1] 37546.92
    print(paste("Anc median, mean: ", median(sim_df_sort_alldat$NANC), ", ", mean(sim_df_sort_alldat$NANC)))    
    #[1] 228071
    #[1] 639442.7
    sim_df_sort_alldat_NE_long <- gather(sim_df_sort_alldat, population, Ne, c(NPOP,NANC), factor_key=TRUE)
    sim_df_sort_alldat_NE_long$Ne <- (sim_df_sort_alldat_NE_long$Ne/2)/100000
    pl <- boxplot(formula = Ne ~ population, data = sim_df_sort_alldat_NE_long, col = cols_fsc[1:2], yaxt = "n",
        names = c(expression(Avg(N[1] + N[2])), expression(N[Ancestral])), border = cols_fsc[1:2], medcol = "gray20",
        ylab = expression("N"[e] ~ (10^5)), ylim = c(0, 7), xlab = NULL)
    axis(2, las = 2)
    #abline(h=22500/100000)
    segments(x0 = 0.5, x1 = 1.5, y0 = (22500/2)/100000, y1=(22500/2)/100000, col = seg_col)
    #abline(h=225000/100000)
    #abline(h=405000/100000)
    segments(x0 = 1.5, x1 = 2.5, y0 = (405000/2)/100000, y1=(405000/2)/100000, col = seg_col)
}
png(file = paste0("/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/panel_sup_fastsimcoal_part_b_ne_sim.png"), 
    res = 800, units = "in", height = 8, width = 7)
par(mar= c(3.2, 4.6, 1, 1))
ne_sim()
dev.off()


##### Split time data dat
split_dat <- function(){
    #print(median(df_sort_alldat$TDIV)); print(mean(df_sort_alldat$TDIV))
    print(paste("splittime median, mean: ", median(df_sort_alldat$TDIV), ", ", mean(df_sort_alldat$TDIV)))    
    #[1] 26893.5
    #[1] 25022.85
    boxplot(df_sort_alldat$TDIV/3/1000, ylab = "Split time (1000 years)", xpd = NA, yaxt = "n",
        col = cols_fsc[3], names = "Population pairs", ylim = c(0,150000/3/1000), medcol = "gray20", border = cols_fsc[3])
    axis(2, las = 2)
    axis(1, at = 1, labels = "Population pairs")
}
png(file = paste0("/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/panel_sup_fastsimcoal_part_c_split_data.png"), 
    res = 800, units = "in", height = 8, width = 5)
par(mar= c(3.2, 4.6, 1, 1))
split_dat()
dev.off()

#dev.off()


##### Split time data SIM
split_sim <- function(){
    #print(median(sim_df_sort_alldat$TDIV)); print(mean(sim_df_sort_alldat$TDIV))
    print(paste("splittime median, mean: ", median(sim_df_sort_alldat$TDIV), ", ", mean(sim_df_sort_alldat$TDIV)))    
    #[1] 39226
    #[1] 40549.71
    boxplot(sim_df_sort_alldat$TDIV/3/1000, ylab = "Split time (1000 years)", xpd = NA, yaxt = "n",
        col = cols_fsc[3], names = "Population pairs", ylim = c(0,150000/3/1000), medcol = "gray20", border = cols_fsc[3])
    axis(2, las = 2)
    axis(1, at = 1, labels = "Population pairs")
    segments(x0 = 0.65, x1 = 1.35, y0 = 45000/3/1000, y1=45000/3/1000, col = seg_col)
    #abline(h=45000)
}

png(file = paste0("/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/panel_sup_fastsimcoal_part_d_split_sim.png"), 
    res = 800, units = "in", height = 8, width = 5)
par(mar= c(3.2, 4.6, 1, 1))
split_sim()
dev.off()


##### Gene flow 
col_pan <- rev(inferno(3, begin = .2, end = .8))
col_pan <- rev(inferno(2, begin = .52, end = .71))
col2rgb(col_pan)
#      [,1] [,2]
#red    194  245
#green   59  123
#blue    79   23
gene_flow_plot <- function(){
    xlimis = c(0,5000);     xpd_is = FALSE;     mgp_is = c(3,1,0) # margin line for ax title, labels and line (3,1,0)
    cex_adj = .9
    
    # migration ancestral
    median_frame <- aggregate(df_sort, list(df_sort$migration_divide), FUN = median, rm.na = T)
    median_frame$mig_anc <- median_frame$mig_anc*100000
    median_frame$mig_rec <- median_frame$mig_rec*100000
    migmax <- max(c(median_frame$mig_rec, median_frame$mig_anc))
    #migmax <- 25
    ylims <- c(min(c(median_frame$mig_rec, median_frame$mig_anc, 0)), migmax)
    #par( new = T)
    plot(median_frame$Group.1, median_frame$mig_anc, #ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), migmax), 
        col = col_pan[1], xlab = NA, frame.plot = F, ylab = expression("Gene flow rate" ~ (10^-5)),
        xlim = xlimis, xpd = xpd_is, type = "l", xaxt = 'n', las = 2,  ylim = ylims)
    axis(1, at = seq(0, 5000, 1000))

    
    # smooth logged and backtransformed
    #func_values <- predict(smooth.spline(x=median_frame$Group.1, y=log(median_frame$mig_rec)), type = "response")
    #func_values <- exp(func_values$y)
    #z_vec <- vector()
    #for(each in seq(NROW(func_values))){
    #    if( each > 2){
    #        last_val <- func_values[each]
    #        sum_f <- sum(func_values[1:(each-1)])
    #        count <- NROW(func_values[1:(each)])
    #        z_val <- last_val * count - sum_f
    #        z_vec <- append(z_vec, z_val)
    #    }else{
    #        z_vec <- append(z_vec, func_values[each])
    #    }
    #}
    #print(NROW(z_vec)); print(NROW(median_frame$Group.1))
    #lines(median_frame$Group.1, (z_vec), col = col_pan[2], lty = 2, ylim = ylims, xlim = xlimis)


    # smooth non-logged
    func_values <- predict(smooth.spline(x=median_frame$Group.1, y=median_frame$mig_rec), type = "response")
    func_values <- func_values$y
    #print(func_values)
    z_vec <- vector()
    count_vec <- vector()
    sumo_vec <- vector()
    last_val_vec <- vector()
    for(each in seq(NROW(func_values))){
        if( each > 2){
            last_val <- func_values[each]
            # last_val = (sum(z_vec) + z_val)/(NROW(z_vec) + 1)
            # last_val = (sum(z_vec) + z_val)/(NROW(z_vec) + 1)
            # last_val * NROW(z_vec + 1) = sum(z_vec ) + z_val
            # last_val * NROW(z_vec + 1) - sum(z_vec ) =  z_val
            sum_o <- sum(z_vec[1:(each-1)])
            count_o <- NROW(z_vec[1:(each)])
            # avg_n <- sum(new_values[1:(each-1)]+last_val)/NROW(new_values[1:(each)])
            z_val <- (last_val * count_o) - sum_o
            z_vec <- append(z_vec, z_val)
            count_vec <- append(count_vec, count_o)
            sumo_vec <- append(sumo_vec, sum_o)
            last_val_vec <- append(last_val_vec, last_val)
        }else{
            z_vec <- append(z_vec, func_values[each])
        }
    }

    print(NROW(z_vec)); print(NROW(median_frame$Group.1))
    lines(median_frame$Group.1, (z_vec), col = col_pan[2], lty = 2, ylim = ylims, xlim = xlimis)
    lines(median_frame$Group.1, func_values, col = col_pan[2], lty = 3, ylim = ylims, xlim = xlimis)


    # regular gaussian logged fit
    #func_values <- predict(glm(median_frame$mig_rec ~ poly(median_frame$Group.1, 2, raw = F), family=gaussian(link = "log")), type = "response")
    #z_vec <- vector()
    #count_vec <- vector()
    #sumo_vec <- vector()
    #last_val_vec <- vector()
    #for(each in seq(NROW(func_values))){
     #   if( each > 2){
      #      last_val <- func_values[each]
       #     # last_val = (sum(z_vec) + z_val)/(NROW(z_vec) + 1)
        #    # last_val = (sum(z_vec) + z_val)/(NROW(z_vec) + 1)
         #   # last_val * NROW(z_vec + 1) = sum(z_vec ) + z_val
          #  # last_val * NROW(z_vec + 1) - sum(z_vec ) =  z_val
           # sum_o <- sum(z_vec[1:(each-1)])
            #count_o <- NROW(z_vec[1:(each)])
#            # avg_n <- sum(new_values[1:(each-1)]+last_val)/NROW(new_values[1:(each)])
 #           z_val <- (last_val * count_o) - sum_o
  #          z_vec <- append(z_vec, z_val)
   #         count_vec <- append(count_vec, count_o)
    #        sumo_vec <- append(sumo_vec, sum_o)
     #       last_val_vec <- append(last_val_vec, last_val)
      #  }else{
       #     z_vec <- append(z_vec, func_values[each])
        #}
    #}
    #print(NROW(z_vec)); print(NROW(median_frame$Group.1))
    #lines(median_frame$Group.1, (z_vec), col = "yellow", lty = 2, ylim = ylims, xlim = xlimis)
    #lines(median_frame$Group.1, func_values, col = "yellow", lty = 3, ylim = ylims, xlim = xlimis)


    # data plotted:
    #mtext(expression("Gene flow" ~ (10^-5)), side = 2, adj = -.7, padj = -3.5, xpd = NA)    #title(ylab = expression("Migration rate" ~ (10^-5)), xpd = NA, mgp = mgp_is, line = 6.5)
    lines(median_frame$Group.1, median_frame$mig_rec, ylim = ylims, 
        col = col_pan[2], xlim = xlimis, xpd = xpd_is)

    title(xlab = "Years bp", xpd = NA, mgp = c(2.5,1,0))
  
    legend("topleft", legend=c( expression("t"[1] ~ "gene flow (Historical)"), expression("t"[2] ~ "gene flow (Recent)"), "Predicted actual gene flow", "Smoothed fit"),
        col = col_pan[c(1,2,2,2)], lty = c(1,1,2,3),  bty = "n", xpd = NA, ncol = 1,
        xjust = 0, text.width = 1000)
}
#png(file = paste0("/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/median_migration_Ne_agri_ONE_FITline_", name_append,"_",dat_rem,".png"), res = 800, units = "in", height = 6, width = 10)
png(file = paste0("/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/panel_sup_fastsimcoal_part_e_geneflow.png"), 
    res = 800, units = "in", height = 6, width = 10)
par(mar = c(4.2, 4.8, 1, 1))
gene_flow_plot()
dev.off()
# gamma glm


#gene_flow_plot()
#dev.off()

## Read in model from pdf

fsc_model_img <- image_read_pdf("/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/Entomobrya_nicoleti/EntNic/figures/fastsimcoal/fastsimvoal_upd.pdf", density = 1000)
#fsc_model_img <- image_crop(fsc_model_img, "1600x1050+1200+650")
fsc_model_img <- image_trim(fsc_model_img)
fsc_model_img <- image_convert(fsc_model_img, matte = TRUE)
#fsc_model_img <- fsc_model_img + theme(plot.margin = margin(-2, -2, -2, -2))
#fsc_model <- rasterGrob(fsc_model_img, interpolate=F)    # Make them into raster objects
#fsc_model <- grid.raster(fsc_model_img, interpolate=F)    # Make them into raster objects
model_img <- function(){
    #plot.new()
    grid.raster(fsc_model_img, x = unit(0.785, "npc"), y = unit(0.18, "npc"), interpolate=F, width = unit(3, "in"), height = unit(2.1, "in"))
}




# plot panel:
# Ne data   -   Ne sim
# Split data    -   Split sim
# gene flow dat     -   model explaining x
    # gene flow without agri and ne, but with mig rec fitted

png(file = paste0("/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/panel_sup_fastsimcoal.png"), 
    res = 800, units = "in", height = 12, width = 8.5)
    par(mar = c(2.5, 5.5, 1, 2), oma = c(2,0,2,0))
    lay_mat <- matrix(   c( 1,1,1,1,2,2,2,2,
                            3,3,3,3,4,4,4,4,
                            5,5,5,5,5,6,6,6),
                            nrow = 3, byrow = T)
    layout(lay_mat)
    ne_data_plot()
    #u <- par("usr")
    #text(u[1], u[4], labels = "a",col = "black", xpd = NA, cex = 1.5, pos = 2)
    mtext("a", side = 3, adj = -.1, cex = 1.2)
    ne_sim()
    mtext("b", side = 3, adj = -.1, cex = 1.2)
    split_dat()
    mtext("c", side = 3, adj = -.1, cex = 1.2)
    split_sim()
    mtext("d", side = 3, adj = -.1, cex = 1.2)
    gene_flow_plot()
    mtext("e", side = 3, adj = -.1, cex = 1.2)
    model_img()
    mtext("f", side = 3, line = -5.5, adj = 1, cex = 1.2)
    mtext("Data SFS", side = 3, outer = T, adj = .26)
    mtext("Simulated SFS", side = 3, outer = T, adj = .83)
dev.off()
