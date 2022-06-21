# BJL
# dataprep unite
# import csv files made by 01_dataprep.R
# these are the pool of each timepoint filtered for coverage
# unites for only bases covered in all timepoints
# calculates percentage methylation 
# saves files

# libraries
library(tidyverse)
library(lubridate)
library(methylKit)

# variables
coverage_threshold <- 5 # set this as appropriate 
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
                        "_pooled_df_filt",
                        coverage_threshold,
                        ".csv"))

tmpnt_pooled_filt <- methRead(location = files,
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


# unite to include only bases covered in all samples
tmpnt_pooled_filt_united <- methylKit::unite(tmpnt_pooled_filt)
#save
saveRDS(tmpnt_pooled_filt_united,
        file = paste0("united_pooled",
                      strand,
                      "_filt",
                      coverage_threshold,
                      ".RDS"))
# make percentage methylation file
pc_tmpnt_pooled_filt_united <- percMethylation(tmpnt_pooled_filt_united,
                                               rowids = T)
#save
saveRDS(pc_tmpnt_pooled_filt_united,
        file = paste0("pc_united_pooled",
                      strand,
                      "_filt",
                      coverage_threshold,
                      ".RDS"))

sessionInfo()
q()



