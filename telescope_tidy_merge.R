# 20210208
# BJL
# merge telescope_report_tidied raw and cpm count files

# packages
library(tidyverse)
library(rtracklayer)

# directories


workingdir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/repeat_expression/"

# raw file
# import files
setwd(workingdir)
files <- list.files(pattern = "tidied.txt")

counts <- files %>%
  set_names() %>% 
  map_dfr(read_delim, delim = "\t", .id = "source")

saveRDS(counts, "telescope_report_tidied_merged.RDS")


# cpm file
# import files
setwd(workingdir)
files <- list.files(pattern = "tidied_cpm.txt")

counts <- files %>%
  set_names() %>% 
  map_dfr(read_delim, delim = "\t", .id = "source")

saveRDS(counts, "telescope_report_tidied_cpm_merged.RDS")