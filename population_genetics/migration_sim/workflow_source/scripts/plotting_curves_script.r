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

# remove pairs with split time below 7000
df_sort <- df_sort[which(df_sort$TDIV > 7000),]
nrow(df_sort)


png(file = "/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/median_migration_curve_EntNic.png", res = 800, units = "in", height = 7, width = 7)
mean_frame <- aggregate(df_sort, list(df_sort$migration_divide), FUN = mean, rm.na = T)
median_frame <- aggregate(df_sort, list(df_sort$migration_divide), FUN = median, rm.na = T)
plot(mean_frame$Group.1, mean_frame$mig_rec, type = "p", ylim = c(min(c(mean_frame$mig_rec, mean_frame$mig_anc)), .0001), 
    ylab = "Migration (ind? chr? per generation)", xlab = "Migration divide (Generations)")
points(mean_frame$Group.1, mean_frame$mig_anc, col = "red", ylim = c(min(c(mean_frame$mig_rec, mean_frame$mig_anc)), .0001))
points(mean_frame$Group.1, median_frame$mig_rec, pch = 2, col = "black", ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), .0001))
points(mean_frame$Group.1, median_frame$mig_anc, pch = 2, col = "red", ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), .0001))
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
points(mean_frame$Group.1, median_frame$mig_rec, pch = 2, col = "black", ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), .00002))
points(mean_frame$Group.1, median_frame$mig_anc, pch = 2, col = "red", ylim = c(min(c(median_frame$mig_rec, median_frame$mig_anc)), .00002))
legend("topright", legend=c("Ancestral migration", "Recent migration"), col = c("red", "black"), pch = 1, bty = "n")
legend("topright", legend=c("Mean", "Median"), col = "black", pch = c(1,2), bty = "n", inset = c(.15,.09))
#plot(mean_frame$Group.1, mean_frame$mig_rec, type = "p", ylim = c(min(c(mean_frame$mig_rec, mean_frame$mig_anc)), max(c(mean_frame$mig_rec, mean_frame$mig_anc))))
#points(mean_frame$Group.1, mean_frame$mig_anc, col = "red", ylim = c(min(c(mean_frame$mig_rec, mean_frame$mig_anc)), max(c(mean_frame$mig_rec, mean_frame$mig_anc))))
dev.off()


par(mfrow=c(1,3))
boxplot(TDIV ~ migration_divide, data = df_sort)
boxplot(mig_anc ~ migration_divide, data = df_sort)
boxplot(mig_rec ~ migration_divide, data = df_sort)
dev.off()



# plot mean demography
pops_to_exclude <- "JEJ-C119|KoeJ-C212|VAJ-C24|HHJ-C162|SBJ-C52|LVJ-C76" #list them, separate by |
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



