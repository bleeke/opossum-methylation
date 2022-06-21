# 20210126
# BJL
# build a repeat annotation
# with useful information from both
# UCSC Table browser '.gtf' output 
# and UCSC Table browser ".txt" output


# packages
library(tidyverse)
library(rtracklayer)

# directories
annotdir <- "/camp/lab/turnerj/working/Bryony/annotations"


# repeats data
#build a repeat annotation with useful information from both UCSC annot formats
setwd(annotdir)

repeats_txt <- read.delim("20210112_UCSC_monDom5_RepeatMasker.txt") # for Class, Family

repeats_gtf <- import("20210123_UCSC_monDom5_RepeatMasker_sorted_shifted_locus.gtf") # for repeat positions and unique loci and corrected Xchr positions.

#  make GRanges; make df with desired attributes; make attribute df into GRanges values
repeats_GR <- GRanges(seqnames = repeats_txt$genoName,
                     ranges = IRanges(start = repeats_txt$genoStart + 1,
                                      end = repeats_txt$genoEnd),
                                      strand = repeats_txt$strand)

tmp_attributes <- repeats_txt %>% select(repName, repClass, repFamily)

values(repeats_GR) <- tmp_attributes
rm(tmp_attributes)

# makes sure order will match
# hard coding the order in seqlevels 
chrom_order <- c(paste0("chr", 1:8), "chrUn", "chrX", "chrM")  
seqlevels(repeats_GR) <- chrom_order
seqlevels(repeats_gtf) <- chrom_order
repeats_GR <- sort(repeats_GR) 
repeats_gtf <- sort(repeats_gtf)


if(isTRUE(all.equal(repeats_GR$repName, as.factor(repeats_gtf$gene_id)))) {
  df_repeats <- data_frame(chrom = as.character(seqnames(repeats_gtf)),
                           start = start(repeats_gtf),
                           end = end(repeats_gtf),
                           strand = as.character(strand(repeats_gtf)),
                           repName = repeats_gtf$gene_id,
                           locus = repeats_gtf$locus,
                           repClass = repeats_GR$repClass,
                           repFamily = repeats_GR$repFamily) 
} else {
  print("objects not equal")
}


# tidy column for repeat class to incorporate "?" repeats into assumed classes
df_repeats$repClass[grepl("DNA?", df_repeats$repClass, fixed = T)] <- "DNA"
df_repeats$repClass[grepl("SINE?", df_repeats$repClass, fixed = T)] <- "SINE"
df_repeats$repClass[grepl("LINE?", df_repeats$repClass, fixed = T)] <- "LINE"
df_repeats$repClass[grepl("LTR?", df_repeats$repClass, fixed = T)] <- "LTR"
df_repeats$repClass[grepl("RC?", df_repeats$repClass, fixed = T)] <- "RC"
df_repeats$repClass[grepl("Unknown?", df_repeats$repClass, fixed = T)] <- "Unknown"


write_delim(df_repeats, "20210126_UCSC_monDom5_repeatmasker_sorted_shifted_classfamily.txt", delim = "\t")