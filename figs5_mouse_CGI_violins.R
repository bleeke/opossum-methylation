# 20220328
# BJL
# violin plots of DNA methylation at CGIs
# mouse allele-specific X chromosome data


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
library(ggridges)
#------------------------------------------------------------------

#------------------------------------------------------------------
# directories
#------------------------------------------------------------------
dataDir <- "/camp/lab/turnerj/working/Bryony/mouse_adult_xci/allele_specific/data/bs-seq/methyl_extract/20190906_methyl_extract"
plotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/plots/"
sampleDir <- "/camp/lab/turnerj/working/Bryony/mouse_adult_xci/allele_specific/data/bs-seq/sample_data"
annotDir <- "/camp/lab/turnerj/working/Bryony/annotations"

# variables
genomes_keep <- c("black6", "spretus")  # genome context of interest; choice of "all" (all mapped reads - #allele_flagged files), "genome1" (mat. files), "genome2", (pat files)

# load data

# sample data
setwd(sampleDir)
sample_data <- read.csv("20190930_sample_data.csv") %>%
  filter(genome %in% genomes_keep) %>%
  mutate(condition = paste0(sex, ".", tissue, ".", genome)) 

# read in previous object generated without male.genome2 to use for X chromosome CpGs (improve #number CpGs kept at "unite")
setwd(dataDir)
X_data <- get(load("2019-12-02fg1_fg2_mg1_allelic_data_merge_locount3_minpergroup2.R"))



# CGIs

setwd(annotDir)
load("2019-11-06mm10_CGI_types_granges.R")
seqlevelsStyle(CGIs) <- "NCBI"

# intersect with methylation file
CGIData <- regionCounts(object = X_data,
                        regions = CGIs,
                        save.db = FALSE) # counts

pcCGIs <- percMethylation(CGIData,
                          rowids = TRUE) # percentages
pcCGIs <- pcCGIs %>% 
  as.data.frame() %>%
  rownames_to_column(var = "position") 


# plots Xa vs Xi vs male X

plotdata <- pcCGIs %>%
  separate(position,
           into = c("chromosome", "start", "end"),
           sep = "[.]",
           remove = FALSE) %>% # adds columns for chromosome, start, end, retains position column
  gather(key = limsid,
         value = pc_methylation,
         -position,
         -chromosome,
         -start,
         -end) %>% #makes data long
  separate(limsid,
           into = c("limsid",
                    "sex",
                    "tissue",
                    "genome"),
           remove = TRUE) %>%
  tidyr::unite(sex, genome,
               col = "condition",
               sep = ".") %>%
  mutate(chr = case_when(chromosome == "X" ~ "chrX",
                         chromosome != "X" ~ "autosome")) %>%
  filter(chr == "chrX")


fills_plot <- c("black", "forestgreen")

plotdata %>%
  ggplot(aes(y = pc_methylation,
           x = condition,
           alpha = factor(condition),
           linetype = factor(condition),
          fill = tissue))+
  geom_violin(scale = "width",
              size = 1.1)+
  theme_bw()+
  xlab("")+
  ylab("CGI methylation")+
  facet_grid(.~tissue)+
  scale_fill_manual(values = fills_plot)+
  scale_linetype_manual(values = c(2, 3, 1))+
  scale_alpha_manual(values = c(0.9, 0.5, 0.1))+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 16),
        legend.position = "none",
        legend.box = "horizontal",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

setwd(plotDir)
ggsave("mouse_allele_specific_Xchr_CGIs.pdf", width = 13, height = 4.5)