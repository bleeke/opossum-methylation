# BJL
# dataprep unite
# import csv files made by 01_dataprep.R
# these are the pool of each timepoint 
# unites to one methylKit object
# keeps all sites 
# saves file

# libraries
library(tidyverse)
library(lubridate)
library(methylKit)

# variables
coverage_threshold <- 1
strand <- "_destranded"

# directories
dataDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/methylKit_out"

# import dfs filtered for coverage using methylkit 
# create list of samples
samples <- as.list(c("1.5_embryo", "2.5_embryo", "3.5_embryo",
                     "4.5_embryo", "5.5_embryo", "6.5_embryo",
                     "7.5_ED", "7.5_embryo", "7.5_TE",
                     "oocyte", "sperm"))

condition <- c(2, 3, 4, 5, 6, 7, 8, 10, 9, 0, 1)

# load data
setwd(dataDir)

files <- as.list(paste0(samples,
                        strand,
                        "_pooled_df",
                        ".csv"))

tmpnt_pooled <- methRead(location = files,
                              sample.id = samples,
                              assembly = "MonDom5",
                              pipeline = list(fraction = TRUE,
                                              chr.col = 1,
                                              start.col = 2,
                                              end.col = 3,
                                              coverage.col = 5,
                                              strand.col = 4,
                                              freqC.col = 8),
                              context = "CpG",
                              resolution = "base",
                              treatment = condition,
                              mincov = coverage_threshold)


# unite to make it a methylBase object but use 'min.per.group=0L' to keep all sites 
tmpnt_pooled_united <- methylKit::unite(tmpnt_pooled,
                                        min.per.group = 0L)
#save
saveRDS(tmpnt_pooled_united,
        file = paste0("united_pooled",
                      strand,
                      ".RDS"))

sessionInfo()
q()



