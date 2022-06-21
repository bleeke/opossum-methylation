# 20220407
# BJL
# allelic X chromosome methylation 
# opossum adult tissues
# summary of genomic features 

# 
#------------------------------------------------------------------
# libraries
#------------------------------------------------------------------
library(methylKit)
library(tidyverse)
library(lubridate)
library(refGenome)
library(GenomicFeatures)
library(genomation)
library(scales)
library(rtracklayer)
#------------------------------------------------------------------

#------------------------------------------------------------------
# directories
#------------------------------------------------------------------
dataDir <- "/camp/lab/turnerj/working/Bryony/opossum_adult/allele-specific/data/bs-seq/methyl_extract"
sampleDir <- "/camp/lab/turnerj/working/Bryony/opossum_adult/allele-specific/data/bs-seq/sample_data"
plotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/plots/"
annotDir <- "/camp/lab/turnerj/working/Bryony/annotations"
rnaDir <- "/camp/lab/turnerj/working/Bryony/opossum_adult/allele-specific/data/rna-seq"
OOPsDir <- "/camp/lab/turnerj/working/shared_projects/OOPs/genome"
#------------------------------------------------------------------

#------------------------------------------------------------------
# variables
#------------------------------------------------------------------
genomes_keep <- c("genome1", "genome2")  # genome context of interest; choice of "all" (all mapped reads - allele_flagged files), "genome1" (mat. files), "genome2", (pat files)
tissues <- c("brain",
             "liver",
             "spleen")
cols_adult <- c("#27408B",
                "#3A5FCD", 
                "#436EEE")
# change index to choose tissue we are working on or use in input script
tissue_keep <- tissues[3]
fill_colour <- cols_adult[3]


#------------------------------------------------------------------

#------------------------------------------------------------------
# load data
#------------------------------------------------------------------

#------------------------------------------------------------------
# sample data 
#------------------------------------------------------------------
# exclude male genome2 as only interested in X chromosome
setwd(sampleDir)
sample_data <- read.csv("20191017_sample_data.csv") %>%
  filter(genome %in% genomes_keep & tissue %in% tissue_keep) %>%
  mutate(condition = paste0(sex, ".", genome)) %>%
  filter(condition != "male.genome2")

#------------------------------------------------------------------
# methylation call files, coverage filtered and merged between samples
#------------------------------------------------------------------
setwd(dataDir)
load(file = paste0(tissue_keep, "allelic_data_merge.R"))

#------------------------------------------------------------------
# annotation objects
#------------------------------------------------------------------

#------------------------------------------------------------------
# CGIs
#------------------------------------------------------------------
setwd(annotDir)
load("2019-11-06mondom5_JZ_mod_shifted_CGI_types_granges.R")
#------------------------------------------------------------------

#------------------------------------------------------------------
# genes
#------------------------------------------------------------------
setwd(OOPsDir)
gtf <- import(file.path(OOPsDir, "Monodelphis_domestica.monDom5.97_mod_shifted_XY_DNMT1.gtf"))
gtf <- gtf[gtf$type == "gene",] # filters for lines type = gene
pseudogenes <- unique(gtf$gene_id[str_detect(gtf$gene_biotype, "pseudogene")])
gtf <- gtf[!gtf$gene_biotype %in% pseudogenes] # exclude pseudogenes
#------------------------------------------------------------------

#------------------------------------------------------------------
# intergenic
#------------------------------------------------------------------
setwd(annotDir)
# bashily - ml SAMtools ; samtools faidx file.fasta ; cut -f1,2 file.faidx > chrom.lengths 
chrom_lengths <- read.delim("mondom5_pseudoY_X-gaps-filled_220819_JZ_chrom_lengths", header = T) # lengths of chromosomes 
chrom_grngs <- GRanges(seqnames = chrom_lengths$chromosome,
                       ranges = IRanges(start = rep(1, times = nrow(chrom_lengths)), end = as.numeric(chrom_lengths$length)))

intergenicRegions <- GenomicRanges::setdiff(chrom_grngs, unstrand(gtf))

rm(chrom_lengths)
#------------------------------------------------------------------

#------------------------------------------------------------------
# promoters
#------------------------------------------------------------------
promoters <- unique(promoters(gtf))
#------------------------------------------------------------------



#------------------------------------------------------------------
# quantify methylation at features
#------------------------------------------------------------------

#------------------------------------------------------------------
# tiles
#------------------------------------------------------------------
tileData <- tileMethylCounts(object = data_merge,
                             win.size = 100,
                             step.size = 100,
                             save.db = FALSE) # counts

pcTiles<- percMethylation(tileData,
                          rowids = TRUE) # percentages
pcTiles <- pcTiles %>% 
  as.data.frame() %>%
  rownames_to_column(var = "position") %>%
  mutate(feature = rep("tile", nrow(pcTiles)))

#------------------------------------------------------------------
# intergenic
#------------------------------------------------------------------
intergenicData <- regionCounts(object = data_merge,
                               regions = intergenicRegions,
                               save.db = FALSE) # counts

pcIntergenic <- percMethylation(intergenicData,
                                rowids = TRUE) # percentages
pcIntergenic <- pcIntergenic %>% 
  as.data.frame() %>%
  rownames_to_column(var = "position") %>%
  mutate(feature = rep("intergenic", nrow(pcIntergenic)))
#------------------------------------------------------------------
# promoters
#------------------------------------------------------------------
promoterData <- regionCounts(object = data_merge,
                             regions = promoters,
                             save.db = FALSE) # counts

pcPromoters <- percMethylation(promoterData,
                               rowids = TRUE) # percentages
pcPromoters <- pcPromoters %>% 
  as.data.frame() %>%
  rownames_to_column(var = "position") %>%
  mutate(feature = rep("promoter", nrow(pcPromoters)))
#------------------------------------------------------------------
# genes
#------------------------------------------------------------------
geneData <- regionCounts(object = data_merge,
                         regions = gtf,
                         save.db = FALSE) # counts
geneData <- geneData[!duplicated(paste0(geneData$chr, geneData$start, geneData$end)),] # remove duplicated gene positions

pcGenes <- percMethylation(geneData,
                           rowids = TRUE) # percentages
pcGenes <- pcGenes %>% 
  as.data.frame() %>%
  rownames_to_column(var = "position") %>%
  mutate(feature = rep("gene", nrow(pcGenes)))

#------------------------------------------------------------------
# CGIs
#------------------------------------------------------------------
CGIData <- regionCounts(object = data_merge,
                        regions = CGIs,
                        save.db = FALSE) # counts

pcCGIs <- percMethylation(CGIData,
                          rowids = TRUE) # percentages
pcCGIs <- pcCGIs %>% 
  as.data.frame() %>%
  rownames_to_column(var = "position") %>%
  mutate(feature = rep("CGI", nrow(pcCGIs)))
#------------------------------------------------------------------
# bind pc methylation data objects together
#------------------------------------------------------------------

pcMeth <- bind_rows(pcTiles,
                    pcIntergenic,
                    pcPromoters,
                    pcGenes,
                    pcCGIs)

#------------------------------------------------------------------
# plot tissue graph
#------------------------------------------------------------------

# format data for graphing
pcMeth <- pcMeth %>% separate(position,
                              into = c("chromosome", "start", "end"),
                              sep = "[.]",
                              remove = FALSE) %>% # adds columns for chromosome, start, end, retains postion column
  mutate(chromosome = str_remove(chromosome, pattern = "chr")) %>%
  filter(chromosome == "X") %>% # filters for chromosome X positions only
  gather(key = limsid,
         value = pc_methylation,
         -position,
         -chromosome,
         -start,
         -end,
         -feature) %>% #makes data long
  separate(limsid,
           into = c("limsid",
                    "sex",
                    "genome"),
           remove = TRUE) %>%
  tidyr::unite(sex, genome,
               col = "condition",
               sep = ".") %>%
  mutate(label = case_when(condition == "male.genome1" ~ "male X",
                           condition == "female.genome1" ~ "female Xa",
                           condition == "female.genome2" ~ "female Xi"))

# setting factor order for plotting
pcMeth$feature <- factor(pcMeth$feature, levels = c("tile",
                                                    "intergenic",
                                                    "CGI", "promoter",
                                                    "gene"))
pcMeth$label <- factor(pcMeth$label, levels = c("male X",
                                                "female Xa",
                                                "female Xi"))


# plot
setwd(plotDir)
pcMeth %>% 
  ggplot(aes(x = feature,
             y = pc_methylation,
             alpha = label,
             linetype = factor(label)))+
  geom_violin(fill = fill_colour,
              scale = "width",
              size = 1.1)+
  scale_linetype_manual(values = c(2, 3, 1))+
  scale_alpha_manual(values = c(0.9, 0.5, 0.1))+
  theme_bw()+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 16),
        legend.position = "none")+
  ylab("methylation level (%%)")+
  xlab("")
ggsave(paste0(tissue_keep, "_allelic_meth_X_opossum.pdf"), height = 3, width = 15)