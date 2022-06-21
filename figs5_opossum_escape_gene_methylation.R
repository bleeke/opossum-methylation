# 20220408
# BJL
# allelic methylation opossum adult tissues
# X chromosome gene bodies
# escape genes



#----------------------------------------------------------------------
# libraries
#----------------------------------------------------------------------
library(methylKit)
library(tidyverse)
library(lubridate)
library(refGenome)
library(GenomicFeatures)
library(genomation)
library(scales)
library(rtracklayer)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# directories 
#----------------------------------------------------------------------
dataDir <- "/camp/lab/turnerj/working/Bryony/opossum_adult/allele-specific/data/bs-seq/methyl_extract"
sampleDir <- "/camp/lab/turnerj/working/Bryony/opossum_adult/allele-specific/data/bs-seq/sample_data"
plotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/plots/"
annotDir <- "/camp/lab/turnerj/working/shared_projects/OOPs/genome/"
escapeDir <- "/camp/home/leekeb/working/Bryony/annotations"
rnaDir <- "/camp/home/leekeb/working/Bryony/opossum_adult/allele-specific/data/rna-seq"
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# variables
#----------------------------------------------------------------------
genomes_keep <- c("genome1", "genome2")  # genome context of interest; choice of "all" (all mapped reads - allele_flagged files), "genome1" (mat. files), "genome2", (pat files)

tissues <- c("brain",
             "liver",
             "spleen")
cols_adult <- c("#27408B",
                "#3A5FCD", 
                "#436EEE")
# change index to choose tissue we are working on or use in input script
tissue_keep <- tissues[2]
fill_colour <- cols_adult[2]
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# import data
#----------------------------------------------------------------------

#------------------------------------------------------------------
# sample data
#------------------------------------------------------------------
setwd(sampleDir)
sample_data <- read.csv("20191017_sample_data.csv") %>%
  filter(genome %in% genomes_keep & tissue %in% tissue_keep) %>%
  mutate(condition = paste0(sex, ".", genome)) 
#------------------------------------------------------------------

#----------------------------------------------------------------------
# methylation call files, coverage filtered and merged between samples
#----------------------------------------------------------------------
setwd(dataDir)
load(file = paste0(tissue_keep, "allelic_data_merge.R"))
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# opossum gtf file
#----------------------------------------------------------------------
setwd(annotDir)
gtf <- import(file.path(annotDir, "Monodelphis_domestica.monDom5.97_mod_shifted_XY.gtf"))
gtf <- gtf[gtf$type == "gene",] # filters for lines type = gene
pseudogenes <- unique(gtf$gene_id[str_detect(gtf$gene_biotype, "pseudogene")])
gtf <- gtf[!gtf$gene_biotype %in% pseudogenes] # exclude pseudogenes
genes <- unique(unstrand(gtf)) # excludes non-unique gene positions # defines unique genes
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# imports expressed genes objects
#----------------------------------------------------------------------
setwd(rnaDir)
# loads genes expressed FPKM > 1 in merged rna-seq
fpkms <- get(load(paste0("merged_", tissue_keep, "_fpkms.R")))
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# defines escape genes
#----------------------------------------------------------------------
# imports BJL defined escape genes for tissue

# imports Wang 2014 defined escape genes
setwd(escapeDir)
wang_escapers <- read.csv("wang_escape_genes.csv")
wang_escapers <- wang_escapers %>% filter(escape_stat == "wang_escaper") # filters out escapers annotated as putative in Wang 2014

escapers  <- wang_escapers %>%
  as_tibble() %>%
  dplyr::select(gene_id, chromosome, start, end)

#----------------------------------------------------------------------

#----------------------------------------------------------------------
# quantify methylation
#----------------------------------------------------------------------
geneData <- regionCounts(object = data_merge,
                         regions = genes,
                         save.db = FALSE) # counts

pcGenes<- percMethylation(geneData,
                          rowids = TRUE) # percentages
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# builds tidy data for plotting
#----------------------------------------------------------------------
plotdat <- as.data.frame(pcGenes) %>%
  rownames_to_column(var = "position") %>%# makes position a columm 
  mutate(position = str_remove(position, pattern = "chr")) %>%
  separate(position,
           into = c("chromosome", "start", "end"),
           sep = "[.]",
           remove = FALSE) %>% # adds columns for chromosome, start, end, retains postion column
  filter(chromosome == "X") %>% # filters for chromosome X positions only
  gather(key = limsid,
         value = pcMeth,
         -position,
         -chromosome,
         -start,
         -end) %>% #makes data long
  separate(limsid,
           into = c("limsid",
                    "sex",
                    "genome"),
           remove = TRUE) %>%
  # removes male
  filter(sex == "female")


# adds column for gene_ids
genes_df <- data.frame(chromosome = seqnames(genes),
                       start = start(genes@ranges),
                       end = end(genes@ranges),
                       gene_id = genes$gene_id) %>%
  mutate(chromosome = str_remove(chromosome,
                                 pattern = "chr")) %>%
  mutate(position = paste0(chromosome,
                           ".",
                           start,
                           ".",
                           end)) %>%
  filter(position %in% plotdat$position) %>%
  arrange(start, end)

#summary(plotdat$position == genes_df$position) # checks row order matches 

plotdat$gene_id <- genes_df$gene_id[which(genes_df$position %in% plotdat$position)]  

plotdat <- plotdat %>%
  filter(gene_id %in% fpkms$gene_id) %>%
  mutate(condition = case_when(
    genome == "genome1" ~ "Xa",
    genome == "genome2" & gene_id %in% escapers$gene_id ~ "Xi\nescaper",
    genome == "genome2" & !gene_id %in% escapers$gene_id ~ "Xi"))

num_regions <- length(unique(plotdat$position))
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# plot escape genes
#----------------------------------------------------------------------
setwd(plotDir)
plotdat %>%
  ggplot(aes(x = condition,
             y = pcMeth,
             alpha = condition,
             linetype = factor(condition)))+
  geom_violin(fill = fill_colour,
              scale = "width",
              size = 1.1)+
  labs(title = paste0("methylation level at ",
                      num_regions,
                      " expressed gene bodies in ",
                      tissue_keep),
       x = NULL,
       y = "methylation (%)",
       name = NULL)+
  scale_linetype_manual(values = c(3, 1, 6))+
  scale_alpha_manual(values = c(0.9, 0.5, 0.1))+
  theme_bw()+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        legend.position = "none")+
  ylim(c(0,100))
ggsave(filename = paste0(tissue_keep, "_allelic_escapers_expressed_genes_violin.pdf"), height = 7, width =5)
#----------------------------------------------------------------------
