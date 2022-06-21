# BJL
# data preparation BS-seq
# import bismark cpg report files for each sample in 'condition' (stage+sex)
# destrands and pools
# writes destranded pooled methylkit object as .RDS file
# extracts a dataframe and filters for minimum coverage = 'coverage_threshold'
# writes a .csv file of destranded pooled coverage filtered data for 'condition' (stage+sex)


# command line arguments
args <- commandArgs(trailingOnly = T)
condition <- args[1]
print(paste0("we are working on condition ", condition))
coverage_threshold <- as.numeric(args[2])
print(paste0("coverage threshold is ", coverage_threshold))


# libraries
library(tidyverse)
library(lubridate)
library(methylKit)
library(rtracklayer)


# directories

dataDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/methyl_extract/bismark_CpG_report"
sampleDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/"
annotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/annotations/"
outDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/methylKit_out/sex_specific/"

# import sample metadata
setwd(sampleDir)
sample_data <- read_delim("PM18153_sample_data.txt", delim = "\t") %>%
  # exclude samples not for analysis
  filter(include_in_analysis == "include")

# import data defined by variable 'condition' from Bismark files
# make metadata for only samples of interest based on condition
condition_samples <- sample_data %>%
  tidyr::unite(col = sample_design_sex, sample_design_char, sample_inferred_sex, sep = "_") %>%
  filter(sample_design_sex == condition)
# define sample names
samples <- as.list(condition_samples$sample_id)
# define files
files <-as.list(paste0(samples,
                       "_merged_fastq.gz_trimmed_bismark_bt2.deduplicated.CpG_report.txt"))
# list the input files for logs
print("we will work on files:")
print(files)
# define sample treatment
design <- condition_samples$sample_design
#import
setwd(dataDir)
dat <- methRead(location = files,
                sample.id = samples,
                assembly = "monDom5",
                pipeline="bismarkCytosineReport",
                treatment = design,
                mincov = 1)


# unite 0L to create object able to be pooled by methylKit
# pool data to create one in silico sample per condition
# unite
dat_united <- methylKit::unite(dat,
                               min.per.group = 0L,
                               destrand = T)
# pool 
dat_pooled <- pool(dat_united,
                   sample.ids = c(condition))
# print nrows for logs
print(paste0("number of rows of destranded united data, i.e. number of CpG sites captured at least once in the pool of ",
             condition,
             " libraries is ",
             nrow(dat_united)))
# save pooled file
setwd(outDir)
saveRDS(dat_pooled,
        file = paste0(condition,
                      "_destranded_pooled.RDS"))
# tidy up
rm(dat, dat_united)



# extract data as df
pooled_data_df <- methylKit::getData(dat_pooled)
pooled_data_df <- pooled_data_df %>% 
  mutate(ratio = numCs1/coverage1)
# save
setwd(outDir)
write_delim(x = pooled_data_df, path = paste0(condition,
                                              "_destranded_pooled_df",
                                              ".csv"), delim = "\t")

# filter by coverage
pooled_data_df_filt <- pooled_data_df %>% filter(coverage1 > (coverage_threshold-1))
#add meth ratio column and save as csv for reading back into methylkit object later
pooled_data_df_filt <- pooled_data_df_filt %>% 
  mutate(ratio = numCs1/coverage1)
# save
setwd(outDir)
write_delim(x = pooled_data_df_filt, path = paste0(condition,
                                                   "_destranded_pooled_df_filt",
                                                   coverage_threshold,
                                                   ".csv"), delim = "\t")
#tidy up
rm(pooled_data_df,
   pooled_data_df_filt)

sessionInfo()
q()
