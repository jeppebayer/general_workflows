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


# plot gns fra stairway likewise
# kun plot cropland
# plot genen flow

# try to combine in one graph?
# setup runs to 14000

png(file = "/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/median_migration_Ne_agri.png", res = 800, units = "in", height = 7, width = 7)
par(mfrow = c(3,1), mar = c(1,5,0,0), oma = c(10, 0, 1, 1))
xlimis = c(0,4000)
xpd_is = FALSE
mgp_is = c(3.5,1,0) # margin line for ax title, labels and line (3,1,0)
# agriculture
par(lwd = 2)

plot(data[which(data$Landuse == col_dat[1,1]),"Year"], data[which(data$Landuse == col_dat[1,1]),"Km2"], type = "l", 
    col = col_dat[1,2], xlab = NA, ylab = NA, xlim = xlimis, xpd = xpd_is, xaxt = 'n', yaxt = 'n')
axis(2, las = 2)
lines(data[which(data$Landuse == col_dat[4,1]),"Year"], data[which(data$Landuse == col_dat[4,1]),"Km2"], col = col_dat[4,2], xlim = xlimis, xpd = xpd_is)
title(ylab = "Percentage of Landcover across Denmark", xlab = "Years bp", xpd = xpd_is)
legend("topright", legend = col_dat[c(1,4),1], col = col_dat[c(1,4),2], lty = 1, bty = "n", xpd = xpd_is)

# migration
mean_frame <- aggregate(df_sort, list(df_sort$migration_divide), FUN = mean, rm.na = T)
median_frame <- aggregate(df_sort, list(df_sort$migration_divide), FUN = median, rm.na = T)
plot(median_frame$Group.1/3, median_frame$mig_anc, ylim = c(min(c(mean_frame$mig_rec, mean_frame$mig_anc)), .00002), 
    col = "red", ylab = NA, xlab = NA,
    xlim = xlimis, xpd = xpd_is, type = "l", xaxt = 'n', yaxt = 'n')
title(ylab = "Migration rate", xpd = xpd_is)
axis(2, las = 2)
lines(median_frame$Group.1/3, median_frame$mig_rec, ylim = c(min(c(mean_frame$mig_rec, mean_frame$mig_anc)), .00002), 
    col = "black", xlim = xlimis, xpd = xpd_is)
legend("topright", legend=c("Ancestral migration", "Recent migration"), col = c("red", "black"), lty = 1, bty = "n", xpd = xpd_is)
#legend("topright", legend=c("Ancestral migration", "Recent migration"), col = c("red", "black"), pch = 1, bty = "n", xpd = xpd_is)
#plot(mean_frame$Group.1, mean_frame$mig_rec, type = "p", ylim = c(min(c(mean_frame$mig_rec, mean_frame$mig_anc)), max(c(mean_frame$mig_rec, mean_frame$mig_anc))))
#points(mean_frame$Group.1, mean_frame$mig_anc, col = "red", ylim = c(min(c(mean_frame$mig_rec, mean_frame$mig_anc)), max(c(mean_frame$mig_rec, mean_frame$mig_anc))))

# Ne trajectories (mean of medians)
ylims = c(min(z_trans_dat), max(z_trans_dat))
year_vec <- ne_dat[,1]*1000
#plot(year_vec,z_trans_dat[,1], type = 'l', col = alpha("orange", .3), ylim = ylims, ylab = "z-transformed Ne", xlab = "Thousand years bp", xlim = xlimis, xpd = xpd_is)
#for(colu in seq(2,ncol(z_trans_dat))){
#    lines(year_vec, z_trans_dat[,colu], col = alpha("orange", .3), ylim = ylims, xlim = xlimis, xpd = xpd_is)
#}
plot(year_vec, rowMeans(z_trans_dat), col = "darkorange", xpd = xpd_is, xlim = xlimis, type = 'l', yaxt = 'n',
    xlab = NA, ylab = NA)
axis(2, las = 2)
title(ylab = "Ne (z-transformed)", xlab = "Years bp", xpd = NA)
#lines(year_vec, rowMeans(z_trans_dat), ylim = ylims, lwd = 5, col = "darkorange", xpd = xpd_is, xlim = xlimis)

dev.off()














for( pop_pair in unique(df_sort$pop_pair)){
    # pop_pair <- "EntNic_aaRJ_C225_vs_EntNic_aeRoe-C36"
    # pop_pair <- "EntNic_aaRJ_C225_vs_EntNic_BJJ-C34"
    # pop_pair <- "EntNic_BIJ-C30_vs_EntNic_BJJ-C34"
    # pop_pair <- "EntNic_BIJ-C30_vs_EntNic_FHJ_C38"
    if( pop_pair == unique(df_sort$pop_pair)[1] ){
        df_sub <- df_sort[which(df_sort$pop_pair == pop_pair),]

        #plot(df_sub$migration_divide, df_sub$mig_rec, ylim = c(min(c(df_sub$mig_rec, df_sub$mig_anc)), .0002))
        ##plot(df_sub$migration_divide, df_sub$mig_rec, ylim = c(min(c(df_sub$mig_rec, df_sub$mig_anc)), max(c(df_sub$mig_rec, df_sub$mig_anc))))
        ##points(df_sub$migration_divide, df_sub$mig_anc, col="red", ylim = c(min(c(df_sub$mig_rec, df_sub$mig_anc)), max(c(df_sub$mig_rec, df_sub$mig_anc))))
        #points(df_sub$migration_divide, df_sub$mig_anc, col="red", ylim = c(min(c(df_sub$mig_rec, df_sub$mig_anc)), .0002))
        #dev.off()

        #plot(df_sub$migration_divide, df_sub$mig_rec, type = "l", ylim = c(min(c(df_sort$mig_rec, df_sort$mig_anc)), max(c(df_sort$mig_rec, df_sort$mig_anc))))
        #lines(df_sub$migration_divide, df_sub$mig_anc, type = "l", col="red", ylim = c(min(c(df_sort$mig_rec, df_sort$mig_anc)), max(c(df_sort$mig_rec, df_sort$mig_anc))))
        
        plot(df_sub$migration_divide, df_sub$mig_rec, type = "l", ylim = c(min(c(df_sort$mig_rec, df_sort$mig_anc)), .0001))
        lines(df_sub$migration_divide, df_sub$mig_anc, type = "l", col="red", ylim = c(min(c(df_sort$mig_rec, df_sort$mig_anc)), .01))

        #plot(df_sub$migration_divide, df_sub$mig_anc, type = "l", ylim = c(min(c(df_sort$mig_rec, df_sort$mig_anc)), .0001))
        #lines(df_sub$migration_divide, df_sub$mig_anc, type = "l", col="red", ylim = c(min(c(df_sort$mig_rec, df_sort$mig_anc)), .01))
    }else{
        df_sub <- df_sort[which(df_sort$pop_pair == pop_pair),]
        #plot(df_sub$migration_divide, df_sub$mig_rec, type = "l", ylim = c(min(df_sort$mig_rec), max(df_sort$mig_rec)))
        #lines(df_sub$migration_divide, df_sub$mig_rec, type = "l", col="black", ylim = c(min(df_sort$mig_rec), max(df_sort$mig_rec)))

        #lines(df_sub$migration_divide, df_sub$mig_rec, type = "l", col="black", ylim = c(min(c(df_sort$mig_rec, df_sort$mig_anc)), max(c(df_sort$mig_rec, df_sort$mig_anc))))
        #lines(df_sub$migration_divide, df_sub$mig_anc, type = "l", col="red", ylim = c(min(c(df_sort$mig_rec, df_sort$mig_anc)), max(c(df_sort$mig_rec, df_sort$mig_anc))))
        
        lines(df_sub$migration_divide, df_sub$mig_rec, type = "l", col="black", ylim = c(min(c(df_sort$mig_rec, df_sort$mig_anc)), .0001))
        lines(df_sub$migration_divide, df_sub$mig_anc, type = "l", col="red", ylim = c(min(c(df_sort$mig_rec, df_sort$mig_anc)), .0001))
        
        #lines(df_sub$migration_divide, df_sub$mig_anc, type = "l", col="red", ylim = c(min(c(df_sort$mig_rec, df_sort$mig_anc)), .0001))
    }
    #dev.off()
}
dev.off()



