#!/usr/bin/env Rscript

# workaround before inputs are ready, and if input is unable to be collected and passed (due to length of argument)

#list_of_files <- list.files(pattern="bestlhoods" , path="/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/2dSFS/2dSFS_EntNic_aaRJ_C225_vs_EntNic_aeRoe-C36_1_vs_2/", recursive = TRUE, full.names = TRUE)
list_of_files <- list.files(pattern="bestlikelyhoods" , path="/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/migration_sim/Collembola/Entomobrya_nicoleti/fsc/2dSFS/", recursive = TRUE, full.names = TRUE)

# remove pops!
# double check, they should already be removed in the pipeline. Better safe than sorry.
list_of_files1 <- list_of_files[grep("JEJ-C119|KoeJ-C212|VAJ-C24|HHJ-C162|SBJ-C52|LVJ-C76|ULJ-C122|SKJ-C16", list_of_files, invert = T)]
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
    pop_pair <- gsub("_[0-9]*_vs_[0-9]*_[0-9]*migDiv_sfs2d.bestlhoods", "", pop_pair)
    
    # make wide format
    data1 <- data.frame(pop_pair, migration_divide, data)
    df <- rbind(df, data1)
}

str(df)
df_sort <- df[order(df[,1], df[,2], method = "auto"),]

# remove pairs with split time below 7000
df_sort <- df_sort[which(df_sort$TDIV > 7000),]



for( pop_pair in unique(df_sort$pop_pair)){
    # pop_pair <- "EntNic_aaRJ_C225_vs_EntNic_aeRoe-C36"
    # pop_pair <- "EntNic_aaRJ_C225_vs_EntNic_BJJ-C34"
    # pop_pair <- "EntNic_GEJ-C42_vs_EntNic_KLJ-C127"
    if( pop_pair == unique(df_sort$pop_pair)[1] ){
        df_sub <- df_sort[which(df_sort$pop_pair == pop_pair),]

        plot(df_sub$migration_divide, df_sub$mig_rec, ylim = c(min(c(df_sub$mig_rec, df_sub$mig_anc)), max(c(df_sub$mig_rec, df_sub$mig_anc))))
        points(df_sub$migration_divide, df_sub$mig_anc, col="red", ylim = c(min(c(df_sub$mig_rec, df_sub$mig_anc)), max(c(df_sub$mig_rec, df_sub$mig_anc))))
        dev.off()

        #plot(df_sub$migration_divide, df_sub$mig_rec, type = "l", ylim = c(min(c(df_sort$mig_rec, df_sort$mig_anc)), max(c(df_sort$mig_rec, df_sort$mig_anc))))
        #lines(df_sub$migration_divide, df_sub$mig_anc, type = "l", col="red", ylim = c(min(c(df_sort$mig_rec, df_sort$mig_anc)), max(c(df_sort$mig_rec, df_sort$mig_anc))))
        
        plot(df_sub$migration_divide, df_sub$mig_rec, type = "l", ylim = c(min(c(df_sort$mig_rec, df_sort$mig_anc)), .01))
        lines(df_sub$migration_divide, df_sub$mig_anc, type = "l", col="red", ylim = c(min(c(df_sort$mig_rec, df_sort$mig_anc)), .01))

        #plot(df_sub$migration_divide, df_sub$mig_anc, type = "l", ylim = c(min(c(df_sort$mig_rec, df_sort$mig_anc)), .0001))
        #lines(df_sub$migration_divide, df_sub$mig_anc, type = "l", col="red", ylim = c(min(c(df_sort$mig_rec, df_sort$mig_anc)), .01))
    }else{
        df_sub <- df_sort[which(df_sort$pop_pair == pop_pair),]
        #plot(df_sub$migration_divide, df_sub$mig_rec, type = "l", ylim = c(min(df_sort$mig_rec), max(df_sort$mig_rec)))
        #lines(df_sub$migration_divide, df_sub$mig_rec, type = "l", col="black", ylim = c(min(df_sort$mig_rec), max(df_sort$mig_rec)))

        #lines(df_sub$migration_divide, df_sub$mig_rec, type = "l", col="black", ylim = c(min(c(df_sort$mig_rec, df_sort$mig_anc)), max(c(df_sort$mig_rec, df_sort$mig_anc))))
        #lines(df_sub$migration_divide, df_sub$mig_anc, type = "l", col="red", ylim = c(min(c(df_sort$mig_rec, df_sort$mig_anc)), max(c(df_sort$mig_rec, df_sort$mig_anc))))
        
        lines(df_sub$migration_divide, df_sub$mig_rec, type = "l", col="black", ylim = c(min(c(df_sort$mig_rec, df_sort$mig_anc)), .01))
        lines(df_sub$migration_divide, df_sub$mig_anc, type = "l", col="red", ylim = c(min(c(df_sort$mig_rec, df_sort$mig_anc)), .01))
        
        #lines(df_sub$migration_divide, df_sub$mig_anc, type = "l", col="red", ylim = c(min(c(df_sort$mig_rec, df_sort$mig_anc)), .0001))
    }
    #dev.off()
}
dev.off()



mean_frame <- aggregate(df_sort, list(df_sort$migration_divide), FUN = mean, rm.na = T)

plot(mean_frame$Group.1, mean_frame$mig_rec, type = "l", ylim = c(min(c(mean_frame$mig_rec, mean_frame$mig_anc)), max(c(mean_frame$mig_rec, mean_frame$mig_anc))))
lines(mean_frame$Group.1, mean_frame$mig_anc, col = "red", ylim = c(min(c(mean_frame$mig_rec, mean_frame$mig_anc)), max(c(mean_frame$mig_rec, mean_frame$mig_anc))))
dev.off()