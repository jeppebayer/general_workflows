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

fst_file <- args_are[1] 
    # fst_file <- "/home/anneaa/EcoGenetics/general_workflow_outputs/population_genetics/fst/Collembola/Entomobrya_nicoleti/EntNic_fst_allPos.fst"
distance_file <- args_are[2]
    # distance_file = "/home/anneaa/EcoGenetics/general_workflows/population_genetics/fst/coordinates.tsv"
outfile_ibd <- args_are[3]
outfile_clad_neig <- args_are[4]
outfile_clad_upgma <- args_are[5]


##############
## Read in data
##############

# fst
fst_dat <- read.table(fst_file, stringsAsFactors = TRUE, na.strings = "na", header = T, sep = '\t')
head(fst_dat)
tail(fst_dat)
dim(fst_dat)
##### OBS CHANGE SEPARATOR FOR POPULATIONS, should be :
# Also remove colnr:
# also bgzip it

fst_mean <- colMeans(as.data.frame(fst_dat), na.rm = TRUE)
head(fst_mean)

# distance
gps_data <- read.table(distance_file, header = TRUE, sep = "\t")
nrow(gps_data)

# subset gps data to tested pops, matching up pop names 
gps_data$"location_code"
if (any(grepl("ae", gps_data$"location_code", fixed = T))){       gps_data$"location_code" <- gsub("ae", "A", gps_data$"location_code")    }
if (any(grepl("aa", gps_data$"location_code",  fixed = T))){       gps_data$"location_code" <- gsub("aa", "A", gps_data$"location_code")    }
if (any(grepl("oe", gps_data$"location_code", fixed = T))){       gps_data$"location_code" <- gsub("oe", "O", gps_data$"location_code")    }
# check if duplicates: Should be true
all( table(gps_data$"location_code") == 1 )

# calculate pairwise distances
vector_km_distances <- vector()
name_pairs_vec <- vector()

for (pop1 in 1:(NROW(pops_letters)-1)){
    gps_row1 <- which(gps_data_sub$"location_code" ==  sapply(strsplit(pops_letters[pop1], "_|-"), "[[", 1))
    for (pop2 in (pop1 + 1):NROW(pops_letters)){
        print(pops_letters[pop1])
        print(pops_letters[pop2])
        name_pairs_vec <- append(name_pairs_vec, paste0(pops_letters[pop1], "-", pops_letters[pop2]))
        gps_row2 <- which(gps_data_sub$"location_code" ==  sapply(strsplit(pops_letters[pop2], "-|_"), "[[", 1))
        distance_km <- distm(c(gps_data_sub$"longitude"[gps_row1], gps_data_sub$"latitude"[gps_row1]), c(gps_data_sub$"longitude"[gps_row2], gps_data_sub$"latitude"[gps_row2]), fun = distHaversine) / 1000
        vector_km_distances <- append(vector_km_distances, distance_km)
        print(distance_km)
    }
}


##############
## Plot Isolation By Distance
##############



