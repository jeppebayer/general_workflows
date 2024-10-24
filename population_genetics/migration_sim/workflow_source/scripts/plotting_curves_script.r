#!/usr/bin/env Rscript

# workaround before inputs are ready, and if input is unable to be collected and passed (due to length of argument)

#list_of_files <- list.files(pattern="bestlhoods" , path="/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/2dSFS/2dSFS_EntNic_aaRJ_C225_vs_EntNic_aeRoe-C36_1_vs_2/", recursive = TRUE, full.names = TRUE)
list_of_files <- list.files(pattern="bestlikelyhoods" , path="/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/2dSFS/", recursive = TRUE, full.names = TRUE)
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
    data <- read.table(file, header = T)
    name <- basename(file)
    count <- NROW(unlist(strsplit(name, "_")))
    migration_divide <- unlist(strsplit(name, "_"))[count-1]
    migration_divide <- as.integer(sapply(strsplit(migration_divide, "mig"), "[[", 1))
    #migration_divide <- migration_divide/3
    # get pop pair
    pop_pair <- gsub("2dSFS_", "", name)
    pop_pair <- gsub("_[0-9]*_vs_[0-9]*_[0-9]*migDiv_sfs2d.bestlikelyhoods", "", pop_pair)
    
    # make wide format
    data1 <- data.frame(pop_pair, migration_divide, data)
    df <- rbind(df, data1)
}

str(df)
df_sort <- df[order(df[,1], df[,2], method = "auto"),]
nrow(df_sort)
# remove pairs with split time below 7000
df_sort <- df_sort[which(df_sort$TDIV > 7000),]
nrow(df_sort)


png(file = "/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/median_migration_curve_EntNic.png", res = 800, units = "in", height = 7, width = 7)
mean_frame <- aggregate(df_sort, list(df_sort$migration_divide), FUN = mean, na.rm = T)
median_frame <- aggregate(df_sort, list(df_sort$migration_divide), FUN = median, rm.na = T)
plot(mean_frame$Group.1, mean_frame$mig_rec, type = "p", ylim = c(min(c(mean_frame$mig_rec, mean_frame$mig_anc)), .0001), 
    ylab = "Migration (ind? chr? per generation)", xlab = "Migration divide (Generations)")
points(mean_frame$Group.1, mean_frame$mig_anc, col = "red", ylim = c(min(c(mean_frame$mig_rec, mean_frame$mig_anc)), .0001))
points(median_frame$Group.1, median_frame$mig_rec, pch = 2, col = "black", ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), .0001))
points(median_frame$Group.1, median_frame$mig_anc, pch = 2, col = "red", ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), .0001))
legend("topright", legend=c("Ancestral migration", "Recent migration"), col = c("red", "black"), pch = 1, bty = "n")
legend("topright", legend=c("Mean", "Median"), col = "black", pch = c(1,2), bty = "n", inset = c(.15,.09))
#plot(mean_frame$Group.1, mean_frame$mig_rec, type = "p", ylim = c(min(c(mean_frame$mig_rec, mean_frame$mig_anc)), max(c(mean_frame$mig_rec, mean_frame$mig_anc))))
#points(mean_frame$Group.1, mean_frame$mig_anc, col = "red", ylim = c(min(c(mean_frame$mig_rec, mean_frame$mig_anc)), max(c(mean_frame$mig_rec, mean_frame$mig_anc))))
dev.off()

png(file = "/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/median_migration_curve_EntNic_zoom.png", res = 800, units = "in", height = 7, width = 7)
mean_frame <- aggregate(df_sort, list(df_sort$migration_divide), FUN = mean, rm.na = T)
median_frame <- aggregate(df_sort, list(df_sort$migration_divide), FUN = median, rm.na = T)
plot(mean_frame$Group.1, mean_frame$mig_rec, type = "p", ylim = c(min(c(mean_frame$mig_rec, mean_frame$mig_anc)), .00002), 
    ylab = "Migration (ind? chr? per generation)", xlab = "Migration divide (Generations)")
points(mean_frame$Group.1, mean_frame$mig_anc, col = "red", ylim = c(min(c(mean_frame$mig_rec, mean_frame$mig_anc)), .00002))
points(median_frame$Group.1, median_frame$mig_rec, pch = 2, col = "black", ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), .00002))
points(median_frame$Group.1, median_frame$mig_anc, pch = 2, col = "red", ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), .00002))
legend("topright", legend=c("Ancestral migration", "Recent migration"), col = c("red", "black"), pch = 1, bty = "n")
legend("topright", legend=c("Mean", "Median"), col = "black", pch = c(1,2), bty = "n", inset = c(.15,.09))
#plot(mean_frame$Group.1, mean_frame$mig_rec, type = "p", ylim = c(min(c(mean_frame$mig_rec, mean_frame$mig_anc)), max(c(mean_frame$mig_rec, mean_frame$mig_anc))))
#points(mean_frame$Group.1, mean_frame$mig_anc, col = "red", ylim = c(min(c(mean_frame$mig_rec, mean_frame$mig_anc)), max(c(mean_frame$mig_rec, mean_frame$mig_anc))))
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
png("/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/Ne_stairway_trajectory_z_and_mean.png", res = 600, height = 7, width = 7, units = "in")
plot(ne_dat[,1], z_trans_dat[,1], type = 'l', col = "grey", ylim = ylims, ylab = "z-transformed Ne (and blue, Nex 10^5)", xlab = "Thousand years bp")
for(colu in seq(2,ncol(z_trans_dat))){
    lines(ne_dat[,1], z_trans_dat[,colu], col = "grey", ylim = ylims)
}
lines(ne_dat[,1], rowMeans(z_trans_dat), ylim = ylims, lwd = 5)
lines(ne_dat[,1], mean_traj, ylim = ylims, lwd = 5, col = "blue")
dev.off()
 # plot only z transformed
ylims = c(min(z_trans_dat), max(z_trans_dat))
png("/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/Ne_stairway_trajectory_z_trans.png", res = 600, height = 7, width = 7, units = "in")
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



# Plot three graphs together
    # Agri, stairway (Ne), migration

col_pan <- rev(plasma(5))
png(file = "/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/median_migration_Ne_agri.png", res = 800, units = "in", height = 7, width = 7)
par(mfrow = c(3,1), mar = c(1.5,5,0,0), oma = c(4, 0, 1, 1))
xlimis = c(0,4000)
xpd_is = FALSE
mgp_is = c(3,1,0) # margin line for ax title, labels and line (3,1,0)
# agriculture
par(lwd = 2)
#range(pretty(c(0, test$y))) 
plot(data[which(data$Landuse == "cropland"),"Year"], data[which(data$Landuse == "cropland"),"Km2"], type = "l", 
    col = col_pan[1], xlab = NA, ylab = NA, xlim = xlimis, xpd = xpd_is, xaxt = 'n', yaxt = 'n', frame.plot = F, ylim = range(pretty(c(0, data[which(data$Landuse == "cropland"),"Km2"]))))
axis(2, las = 2)
#axis(1, labels = F, col = "grey")
lines(data[which(data$Landuse == "urban_area"),"Year"], data[which(data$Landuse == "urban_area"),"Km2"], col = col_pan[2], xlim = xlimis, xpd = xpd_is)
title(ylab = "Percentage Landcover", xpd = xpd_is, mgp = mgp_is)
legend("right", legend = c("cropland", "urban_area"), col = col_pan[1:2], lty = 1, bty = "n", xpd = xpd_is)

# migration
median_frame <- aggregate(df_sort, list(df_sort$migration_divide), FUN = median, rm.na = T)
median_frame$mig_anc <- median_frame$mig_anc*100000
median_frame$mig_rec <- median_frame$mig_rec*100000
#migmax <- 
migmax <- max(c(median_frame$mig_rec, median_frame$mig_anc))
plot(median_frame$Group.1/3, median_frame$mig_anc, #ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), migmax), 
#plot(median_frame$Group.1/3, median_frame$mig_anc, ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), .00002), 
    col = col_pan[3], ylab = NA, xlab = NA, frame.plot = F, 
    xlim = xlimis, xpd = xpd_is, type = "l", xaxt = 'n', yaxt = 'n', ylim = range(pretty(c(0, migmax))) )
title(ylab = expression("Migration rate" ~ (10^-5)), xpd = xpd_is, mgp = mgp_is)
axis(2, las = 2)
#axis(1, labels = F, col = "grey")
lines(median_frame$Group.1/3, median_frame$mig_rec, ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), migmax), 
#lines(median_frame$Group.1/3, median_frame$mig_rec, ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), .00002), 
    col = col_pan[4], xlim = xlimis, xpd = xpd_is)
legend("right", legend=c("Ancestral ", "Recent "), col = col_pan[3:4], lty = 1, bty = "n", xpd = xpd_is)

# Ne trajectories (mean of medians)
year_vec <- ne_dat[,1]*1000
year_vec <- year_vec[which(year_vec < 4500)]
ylims = range(pretty(c(min(rowMeans(z_trans_dat[1:NROW(year_vec),])), max(rowMeans(z_trans_dat[1:NROW(year_vec),])))))
#ylims = range(pretty(c(min(rowMeans(z_trans_dat))), max(rowMeans(z_trans_dat))))


plot(year_vec, rowMeans(z_trans_dat[1:NROW(year_vec),]), col = col_pan[5], xpd = xpd_is, xlim = xlimis, type = 'l', yaxt = 'n',
    xlab = NA, ylab = NA, frame.plot = F, ylim = ylims)
axis(2, las = 2)
#axis(2, las = 2, at = c(-.5, 0, 0.5, 1), xpd = NA)
title(ylab = "Ne (z-transformed)", xpd = NA, mgp = mgp_is)
title(xlab = "Years bp", xpd = NA, mgp = c(3,1,0))

segments( 8000/3, -.8, 8000/3, 4.1, lty = 2, col = "grey", xpd = NA, lwd = 1)
segments( 600/3, -.8, 600/3, 4.1, lty = 2, col = "grey", xpd = NA, lwd = 1)
segments( 100/3, -.8, 100/3, 4.1, lty = 2, col = "grey", xpd = NA, lwd = 1)
dev.off()








# try to combine in one graph?
    # Agri, stairway (Ne), migration
col_pan <- rev(inferno(3, begin = .2, end = .8))
png(file = "/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/median_migration_Ne_agri_ONE.png", res = 800, units = "in", height = 6, width = 10)
par(mfrow = c(1,1), mar = c(1.5,5,0,0), oma = c(6, 4, 1, 1))
xlimis = c(0,4000);     xpd_is = FALSE;     mgp_is = c(3,1,0) # margin line for ax title, labels and line (3,1,0)
cex_adj = .9
# agriculture
par(lwd = 2)
plot(data[which(data$Landuse == "cropland"),"Year"], data[which(data$Landuse == "cropland"),"Km2"], type = "l", 
    col = "white", xlab = NA, ylab = NA, xlim = xlimis, xpd = xpd_is, xaxt = 'n', yaxt = 'n', frame.plot = F, 
    ylim = range(pretty(c(0, data[which(data$Landuse == "cropland"),"Km2"]))))
axis(2, las = 2, line = 0, col = col_pan[1], cex.axis = cex_adj)
rect( 7000/3, -10, 9000/3, 80, col = "grey90", border = NULL, lwd = 0, xpd = F)
rect( 100, -10, 600, 80, col = "grey90", border = NULL, lwd = 0, xpd = F)
lines(data[which(data$Landuse == "cropland"),"Year"], data[which(data$Landuse == "cropland"),"Km2"], col = col_pan[1], xlim = xlimis, xpd = xpd_is)

# migration
median_frame <- aggregate(df_sort, list(df_sort$migration_divide), FUN = median, rm.na = T)
median_frame$mig_anc <- median_frame$mig_anc*100000
median_frame$mig_rec <- median_frame$mig_rec*100000
migmax <- max(c(median_frame$mig_rec, median_frame$mig_anc))
par( new = T)
plot(median_frame$Group.1/3, median_frame$mig_anc, #ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), migmax), 
    col = col_pan[2], ylab = NA, xlab = NA, frame.plot = F, 
    xlim = xlimis, xpd = xpd_is, type = "l", xaxt = 'n', yaxt = 'n', ylim = range(pretty(c(0, migmax))) )
axis(2, las = 2, line = 3, col = col_pan[2], cex.axis = cex_adj)
#mtext(expression("Gene flow" ~ (10^-5)), side = 2, adj = -.7, padj = -3.5, xpd = NA)    #title(ylab = expression("Migration rate" ~ (10^-5)), xpd = NA, mgp = mgp_is, line = 6.5)


lines(median_frame$Group.1/3, median_frame$mig_rec, ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), migmax), 
    col = col_pan[2], xlim = xlimis, xpd = xpd_is)
points(median_frame$Group.1/3, median_frame$mig_rec, ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), migmax), 
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

text(-330, -.78, labels=c("Percentage\n Landcover"), col = col_pan[1], xpd = NA, srt = 90, adj = 0)
text(-600, -.78, labels=c(expression("Gene flow" ~ (10^-5))), col = col_pan[2], xpd = NA, srt = 90, adj = 0)
text(-900, -.78, labels=c("Ne\n(z-transformed)"), col = col_pan[3], xpd = NA, srt = 90, adj = 0)

legend("bottom", legend=c( "Ancestral gene flow", "Recent gene flow", "Cropland cover", "Ne trajectory"),
    col = col_pan[c(2,2,1,3)], lty = 1, pch = c(NA, 19, NA, NA), bty = "n", xpd = NA, inset = -.33, ncol = 2)
dev.off()









# trying to dissolve flat migration line

png(file = "/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/median_migration_Ne_recent_zoom.png", res = 800, units = "in", height = 7, width = 7)
par(mfrow = c(4,1), mar = c(1,5,0,0), oma = c(4, 0, 1, 1))
xlimis = c(0,4000)
xpd_is = FALSE
mgp_is = c(4,1,0) # margin line for ax title, labels and line (3,1,0)
# agriculture
par(lwd = 2)
# migration log
median_frame <- aggregate(df_sort, list(df_sort$migration_divide), FUN = median, rm.na = T)
median_frame$mig_anc <- log(median_frame$mig_anc)
median_frame$mig_rec <- log(median_frame$mig_rec)
migmax <- max(c(median_frame$mig_rec, median_frame$mig_anc))
plot(median_frame$Group.1/3, median_frame$mig_anc, ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), migmax), 
    col = "red", ylab = NA, xlab = NA,
    xlim = xlimis, xpd = xpd_is, type = "l", xaxt = 'n', yaxt = 'n')
title(ylab = "log Migration rate", xpd = xpd_is, mgp = mgp_is)
axis(2, las = 2)
lines(median_frame$Group.1/3, median_frame$mig_rec, ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), migmax), 
    col = "black", xlim = xlimis, xpd = xpd_is)
legend("topright", legend=c("Ancestral migration", "Recent migration"), col = c("red", "black"), lty = 1, bty = "n", xpd = xpd_is)

# migration norm
median_frame <- aggregate(df_sort, list(df_sort$migration_divide), FUN = median, rm.na = T)
median_frame$mig_anc <- median_frame$mig_anc*100000
median_frame$mig_rec <- median_frame$mig_rec*100000
migmax <- max(c(median_frame$mig_rec, median_frame$mig_anc))
plot(median_frame$Group.1/3, median_frame$mig_anc, ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), migmax), 
    col = "red", ylab = NA, xlab = NA,
    xlim = xlimis, xpd = xpd_is, type = "l", xaxt = 'n', yaxt = 'n')
title(ylab = expression("Migration rate" ~ (10^-5)), xpd = xpd_is, mgp = mgp_is)
axis(2, las = 2)
lines(median_frame$Group.1/3, median_frame$mig_rec, ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), migmax), 
    col = "black", xlim = xlimis, xpd = xpd_is)
legend("topright", legend=c("Ancestral migration", "Recent migration"), col = c("red", "black"), lty = 1, bty = "n", xpd = xpd_is)


# Migration recent
median_frame <- aggregate(df_sort, list(df_sort$migration_divide), FUN = median, rm.na = T)
median_frame$mig_anc <- median_frame$mig_anc*100000
median_frame$mig_rec <- median_frame$mig_rec*100000
migmax <- .001
migmin <- .00005
plot(median_frame$Group.1/3, median_frame$mig_rec, ylim = c(migmin, migmax), 
    col = "grey", ylab = NA, xlab = NA,
    xlim = xlimis, xpd = xpd_is, type = "l", xaxt = 'n', yaxt = 'n')
title(ylab = expression("Migration rate" ~ (10^-5)), xpd = xpd_is, mgp = mgp_is)
axis(2, las = 2)
legend("topright", legend=c( "Recent migration"), col = c("grey"), lty = 1, bty = "n", xpd = xpd_is)

# Migration recent
median_frame <- aggregate(df_sort, list(df_sort$migration_divide), FUN = median, rm.na = T)
median_frame$mig_anc <- median_frame$mig_anc*100000
median_frame$mig_rec <- median_frame$mig_rec*100000
migmax <- .00001
plot(median_frame$Group.1/3, median_frame$mig_rec, ylim = c(min(median_frame$mig_rec), migmax), 
    col = "grey", ylab = NA, xlab = NA,
    xlim = xlimis, xpd = xpd_is, type = "l", xaxt = 'n', yaxt = 'n')
title(ylab = expression("Migration rate" ~ (10^-5)), xpd = xpd_is, mgp = mgp_is)
axis(2, las = 2)
legend("topright", legend=c( "Recent migration"), col = c("grey"), lty = 1, bty = "n", xpd = xpd_is)

dev.off()
