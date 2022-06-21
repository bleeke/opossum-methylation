# BJL
# get overlaps of genomic features with CpGs captured at each timepoint of adult opossum BSseq libraries

# libraries
library(tidyverse)
library(methylKit)
library(GenomicRanges)
library(genomation)
library(rtracklayer)


# command line arguments
args <- commandArgs(trailingOnly = T)
condition <- args[1]
print(paste0("we are working on condition ", condition))

# directories
dataDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/adult_bsseq/methylKit_out"
annotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/annotations/"
plotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/plots/"


# get annotation objects
setwd(annotDir)

# repeats
repeats <- read_delim("monDom5_UCSC_RepeatMasker_Xshifted_ClassFamily.txt", delim = "\t", col_names = T)

repeats_gr <- GRanges(seqnames = repeats$chrom,  # coerce to GRanges
                      ranges = IRanges(start = repeats$start, end = repeats$end),
                      class = repeats$repClass)

repeats <- repeats_gr
rm(repeats_gr)


# CGIs
CGIs <- readRDS("mondom5_CGI_Xshifted.RDS")

# genes  
genes <- import("mondom5.97_Xshifted_Y_DNMT1.gtf")
genes <- genes[genes$type == "gene",] # filters for lines type = gene


promoters <- unique(promoters(genes))

# chromosomes
# chrom lengths information extracted via SAMtools ; samtools faidx file.fasta ; cut -f1,2 file.faidx 
chrom_lengths <- read.delim("genome/chrom_lengths", header = T) 
chromosomes <-  GRanges(seqnames = chrom_lengths$chromosome,
                        ranges = IRanges(start = rep(1, times = nrow(chrom_lengths)),
                                         end = as.numeric(chrom_lengths$length)),
                        chromosome = chrom_lengths$chromosome)

# get data object
setwd(dataDir)
data <-read_delim(paste0(condition,  "_destranded_pooled_df_filt5.csv"), delim = "\t")
# make GRanges
data <- GRanges(seqnames = data$chr,
                ranges = IRanges(start = data$start,
                                 end = data$end),
                strand = data$strand)


# get overlaps between CpGs and annots

#repeats
repeat_overlaps <- findOverlaps(repeats, data)

no_repeat <- length(data) - length(repeat_overlaps)
Simple_repeat <- length(repeat_overlaps[repeats$class[queryHits(repeat_overlaps)] == "Simple_repeat",])
Low_complexity <- length(repeat_overlaps[repeats$class[queryHits(repeat_overlaps)] == "Low_complexity",])
LINE <- length(repeat_overlaps[repeats$class[queryHits(repeat_overlaps)] == "LINE",])
SINE <- length(repeat_overlaps[repeats$class[queryHits(repeat_overlaps)] == "SINE",])
DNA <- length(repeat_overlaps[repeats$class[queryHits(repeat_overlaps)] == "DNA",])
LTR <- length(repeat_overlaps[repeats$class[queryHits(repeat_overlaps)] == "LTR",])
scRNA <- length(repeat_overlaps[repeats$class[queryHits(repeat_overlaps)] == "scRNA",])
tRNA <- length(repeat_overlaps[repeats$class[queryHits(repeat_overlaps)] == "tRNA",])
snRNA <- length(repeat_overlaps[repeats$class[queryHits(repeat_overlaps)] == "snRNA",])
Unknown <- length(repeat_overlaps[repeats$class[queryHits(repeat_overlaps)] == "Unknown",])
rRNA <- length(repeat_overlaps[repeats$class[queryHits(repeat_overlaps)] == "rRNA",])
RC <- length(repeat_overlaps[repeats$class[queryHits(repeat_overlaps)] == "RC",])
Satellite <- length(repeat_overlaps[repeats$class[queryHits(repeat_overlaps)] == "Satellite",])

repeat_df <- data.frame(Simple_repeat, Low_complexity,
                        LINE, SINE, DNA, LTR, scRNA, 
                        tRNA, snRNA, Unknown, rRNA,
                        RC, Satellite, no_repeat) %>% 
  gather(key = type, value = count) %>%
  mutate(perc = count/sum(count)*100)

# CGIs
CGI_overlaps <- findOverlaps(CGIs, data)

no_CGI <- length(data) - length(repeat_overlaps)
prom_CGI <- length(CGI_overlaps[CGIs$CGI_type[queryHits(CGI_overlaps)] == "promoter",])
intergenic_CGI <- length(CGI_overlaps[CGIs$CGI_type[queryHits(CGI_overlaps)] == "intergenic",])
intragenic_CGI <- length(CGI_overlaps[CGIs$CGI_type[queryHits(CGI_overlaps)] == "intragenic",])
CGI_df <- data.frame(prom_CGI, intergenic_CGI, intragenic_CGI, no_CGI) %>%
  gather(key = type, value = count) %>%
  mutate(perc = count/sum(count)*100)

# genes and promoters

promoter_overlaps <- findOverlaps(promoters, data) # CpGs that overlap promoters
data_noproms <- subsetByOverlaps(data, promoters, invert = T) # remove the CpGs that overlap promoters before checking how many overlap genes
gene_overlaps <- findOverlaps(genes, data_noproms)

intergenic <- length(data) - (length(promoter_overlaps) + length(gene_overlaps))
promoter <- length(promoter_overlaps)
intragenic <- length(gene_overlaps)
genes_df <- data.frame(intergenic, intragenic, promoter) %>%
  gather(key = type, value = count) %>%
  mutate(perc = count/sum(count)*100)

# chromosomes

chrom_overlaps <- findOverlaps(chromosomes, data)
chr1 <- length(chrom_overlaps[chromosomes$chromosome[queryHits(chrom_overlaps)] == "chr1",])
chr2 <- length(chrom_overlaps[chromosomes$chromosome[queryHits(chrom_overlaps)] == "chr2",])
chr3 <- length(chrom_overlaps[chromosomes$chromosome[queryHits(chrom_overlaps)] == "chr3",])
chr4 <- length(chrom_overlaps[chromosomes$chromosome[queryHits(chrom_overlaps)] == "chr4",])
chr5 <- length(chrom_overlaps[chromosomes$chromosome[queryHits(chrom_overlaps)] == "chr5",])
chr6 <- length(chrom_overlaps[chromosomes$chromosome[queryHits(chrom_overlaps)] == "chr6",])
chr7 <- length(chrom_overlaps[chromosomes$chromosome[queryHits(chrom_overlaps)] == "chr7",])
chr8 <- length(chrom_overlaps[chromosomes$chromosome[queryHits(chrom_overlaps)] == "chr8",])
chrMT <- length(chrom_overlaps[chromosomes$chromosome[queryHits(chrom_overlaps)] == "chrMT",])
chrX <- length(chrom_overlaps[chromosomes$chromosome[queryHits(chrom_overlaps)] == "chrX",])
chrUn <- length(chrom_overlaps[chromosomes$chromosome[queryHits(chrom_overlaps)] == "chrUn",])
chr_df <- data.frame(chr1, 
                     chr2, 
                     chr3, 
                     chr4, 
                     chr5, 
                     chr6, 
                     chr7, 
                     chr8,
                     chrMT,
                     chrX, 
                     chrUn) %>%
  gather(key = type, value = count) %>%
  mutate(perc = count/sum(count)*100)

# bind dfs
df <- bind_rows(list(genes_df = genes_df,
                     chr_df = chr_df,
                     CGI_df = CGI_df,
                     repeat_df = repeat_df), .id = "context")

# save summary df
setwd(dataDir)
saveRDS(df, file = paste0(condition, "_genomic_overlaps.RDS"))
