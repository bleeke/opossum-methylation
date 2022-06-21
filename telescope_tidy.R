# 20210208
# BJL
# tidies telescope output
# calculates cpm
# writes output file for raw and cpm

# packages
library(tidyverse)


# arguments as defined in the bash submission script
args <- commandArgs(trailingOnly = T)
input <- args[1]
workingDir <- args[2]

setwd(workingDir)

# read input file header and reformat to dataframe
runinfo <- read_lines(input, n_max = 1)
runinfo <- data_frame(info = unlist(str_split(runinfo, pattern = "\t")))
runinfo <- runinfo[2:nrow(runinfo),] %>% 
  separate(info, sep = ":", into = c("category", "value")) %>%
  spread(key = "category", value = "value") %>%
  mutate(total_mapped = as.integer(total_fragments)-as.integer(unmapped)) # calculate total mapped reads

# read input file data and calculate cpm
repeatCounts <- read_delim(input, col_names = TRUE, skip = 1, delim = "\t") %>%
  mutate(cpm = (final_count/runinfo$total_mapped)*1000000)

# write to file

# raw

output <- repeatCounts %>%
  select(transcript, final_count)

write_delim(output, paste0(input, "_tidied.txt"), delim = "\t")


# cpm

output <- repeatCounts %>%
  select(transcript, cpm)

write_delim(output, paste0(input, "_tidied_cpm.txt"), delim = "\t")