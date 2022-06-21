# BJL
# mouse allele-specific BSseq


#------------------------------------------------------------------
# libraries
#------------------------------------------------------------------
library(methylKit)
library(tidyverse)
library(lubridate)
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

# data preparation
# done once

#-------
## variables
#genomes_keep <- c("black6", "spretus")  # genome context of interest; choice of "all" (all #mapped reads - #allele_flagged files), "genome1" (mat. files), "genome2", (pat files)
## load data
## sample data
#setwd(sampleDir)
#sample_data <- read.csv("20190930_sample_data.csv") %>%
#  filter(genome %in% genomes_keep) %>%
#  mutate(condition = paste0(sex, ".", tissue, ".", genome)) 
#
## read in previous object generated without male.genome2 to use for X chromosome CpGs (improve #number CpGs kept at "unite")
#setwd(dataDir)
#X_data <- get(load("2019-12-02fg1_fg2_mg1_allelic_data_merge_locount3_minpergroup2.R"))
#
#pcXData <- percMethylation(X_data,
#                           rowids = TRUE)
#pcXData <- pcXData %>% 
#  as.data.frame() %>%
#  rownames_to_column(var = "position") %>%
#  separate(position,
#           into = c("chromosome", "start", "end"),
#           sep = "[.]",
#           remove = FALSE) %>%
#  mutate(chromosome = str_remove(chromosome, pattern = "chr")) %>%
#  filter(chromosome == "X") 

##read in methylation object that includes male.genome2 to use for autosomes
#auto_data <- get(load("2019-10-07allelic_data_merge_locount3_minpergroup2.R")) 
#auto_data <- reorganize(auto_data,
#                        sample.ids = auto_data@sample.ids,
#                        treatment = as.numeric(as.factor(sample_data$condition)),
#                        save.db = FALSE)
#
#pcAutoData <- percMethylation(auto_data,
#                           rowids = TRUE)
#pcAutoData <- pcAutoData %>% 
#  as.data.frame() %>%
#  rownames_to_column(var = "position") %>%
#  separate(position,
#           into = c("chromosome", "start", "end"),
#           sep = "[.]",
#           remove = FALSE) %>%
#  mutate(chromosome = str_remove(chromosome, pattern = "chr")) %>%
#  filter(chromosome != "X")
#
## bind autosomes and X
#plotdata <- bind_rows(pcXData, pcAutoData)
#rm(auto_data, X_data)
#saveRDS(plotdata, "20220327_mouse_allelic_data_merge_locount3_minpergroup2_Xand_auto.RDS")
#-------

#setwd(dataDir)
#plotdata <- readRDS("20220327_mouse_allelic_data_merge_locount3_minpergroup2_Xand_auto.RDS")

# format data for graphing
#pcData <- plotdata %>% separate(position,
#                              into = c("chromosome", "start", "end"),
#                              sep = "[.]",
#                              remove = FALSE) %>% # adds columns for chromosome, start, end, retains position column
#  mutate(chromosome = str_remove(chromosome, pattern = "chr")) %>%
#  gather(key = limsid,
#         value = pc_methylation,
#         -position,
#         -chromosome,
#         -start,
#         -end) %>% #makes data long
#  separate(limsid,
#           into = c("limsid",
#                    "sex",
#                    "tissue",
#                    "genome"),
#           remove = TRUE) %>%
#  tidyr::unite(sex, genome,
#               col = "condition",
#               sep = ".") %>%
#  mutate(chr = case_when(chromosome == "X" ~ "chrX",
#                                chromosome != "X" ~ "autosome")) %>%
#  mutate(chr = factor(chr,
#                      levels = c("autosome",
#                                 "chrX")))
setwd(dataDir)
#saveRDS(pcData, "20220327_mouse_allelic_data_merge_locount3_minpergroup2_Xand_auto_plotreformatted.RDS")

pcData <- readRDS("20220327_mouse_allelic_data_merge_locount3_minpergroup2_Xand_auto_plotreformatted.RDS")




# plot
setwd(plotDir)

fills_plot <- c("black", "forestgreen")


ridge_plot <- pcData %>%
  ggplot(aes(x = pc_methylation,
             y = condition,
             colour = tissue,
             fill = tissue,
             group = condition, 
             linetype = condition)) +
  geom_density_ridges(alpha = 0.2,
                      size = 1.2,
                      show.legend = F)+
  theme_bw()+
  xlab("CpG methylation")+
  ylab("density")+
  facet_grid(chr ~ tissue)+
  scale_fill_manual(values = fills_plot)+
  scale_colour_manual(values = fills_plot)+
  scale_linetype_manual(values = c(2, 1, 3, 4))+
  scale_x_continuous(breaks = c(0, 50, 100),
                     labels = c("0", "0.5", "1"))+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 16),
        legend.position = "top",
        legend.box = "horizontal",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key.size = unit(4,"line"))


setwd(plotDir)
ggsave("ridgeplot_mouse_adult.pdf", plot = ridge_plot, width = 13, height = 5)


# count number of CpGs for legend
count_CpG <- pcData %>% 
  group_by(condition, tissue, chromosome) %>%
  summarise(count = n())
setwd(dataDir)
write.table(count_CpG, "Count_CpG_mouse_allelic_ridges.txt")

