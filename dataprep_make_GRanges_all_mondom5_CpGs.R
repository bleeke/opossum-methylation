# BJL
# make Granges object of every CpG in the opossum genome

# libraries
library(Biostrings)
library(tidyverse)
library(GenomicRanges)

# directories
dataDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/methylKit_out"
annotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/annotations/"

# get opossum genome sequence
setwd(annotDir)
genome_seq <- readDNAStringSet("genome/bryonys_dream.fasta")
# fix chromosomes names
tmp <- names(genome_seq) %>%
  str_remove(pattern = " .*")
names(genome_seq) <- tmp
rm(tmp)


# get all CpGs as a GRanges object, thanks to https://support.bioconductor.org/p/95239/
chrs <- names(genome_seq)[1:12]
cgs <- lapply(chrs, function(x) start(matchPattern("CG", genome_seq[[x]])))
cpgr <- do.call(c, lapply(1:12, function(x) GRanges(names(genome_seq)[x], IRanges(cgs[[x]], width = 2))))

# total number of CpG sites (destranded)
length(cpgr)

# save
setwd(annotDir)
saveRDS(cpgr, file = "all_mondom_cpgs.RDS")