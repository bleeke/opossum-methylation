library(tidyverse)

setwd("/camp/lab/turnerj/working/Bryony/annotations/")

# this is a tidyverse function (preserves quotes where they are needed):
gtf <- read_delim("20210112_UCSC_monDom5_RepeatMasker_sorted.gtf", delim = "\t", col_names = F)
colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")


# 177 repeats are affected:
gtf %>%
	filter(seqname == "chrX", start > 35623969, end < 35710776) %>% 
	count()


# identify repeat rows to shift depending on coordinates:
# 2 repeats btw gaps 1 and 2:
set1 <- c((gtf$seqname == "chrX") & (gtf$start > 35623969) & (gtf$end < 35636737))
# 2 repeats btw gaps 2 and 3:
set2 <- c((gtf$seqname == "chrX") & (gtf$start > 35637245) & (gtf$end < 35638911))
# 120 repeats btw gaps 3 and 4:
set3 <- c((gtf$seqname == "chrX") & (gtf$start > 35650054) & (gtf$end < 35690788))
# 53 repeats btw gaps 4 and 5:
set4 <- c((gtf$seqname == "chrX") & (gtf$start > 35691812) & (gtf$end < 35710776))
# repeats downstream of gap 5:
set5 <- c((gtf$seqname == "chrX") & (gtf$start > 35711969))

# shift repeats according to coordinates:
gtf[set1, ][c("start", "end")] <- gtf[set1, ][c("start", "end")] - 1040
gtf[set2, ][c("start", "end")] <- gtf[set2, ][c("start", "end")] - 1460
gtf[set3, ][c("start", "end")] <- gtf[set3, ][c("start", "end")] + 1207
gtf[set4, ][c("start", "end")] <- gtf[set4, ][c("start", "end")] + 1195
gtf[set5, ][c("start", "end")] <- gtf[set5, ][c("start", "end")] + 1200

# reformat numbers so that they don't appear in scientific notation:
gtf_final <- mutate_if(gtf, is.numeric, as.integer)

# write new gtf to file:
write.table(gtf_final, file = "20210112_UCSC_monDom5_RepeatMasker_sorted_shifted.gtf", sep = "\t", quote = F, col.names = F, row.names = F)
