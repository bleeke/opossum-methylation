# # figure 3 and s4
# diffMeth regions 
# stacked bar plots
# density plots

# libraries
library(tidyverse)
library(plyr)
library(methylKit)
library(rtracklayer)
library(genomation)
library(ggridges)


# directories

dataDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/methylKit_out/"
adultDir  <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/adult_bsseq/methylKit_out/"
annotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/annotations/"
plotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/plots"


# load annotations 
setwd(annotDir)
# load from gtf to GRanges
genes <- import("mondom5.97_Xshifted_Y_DNMT1.gtf")
# filter out entries for transcripts etc and include pseudoY chromosome CDS 'genes'
genes <- genes[genes$type == "gene" | seqnames(genes) == "chrY"]
# exclude pseudogenes
pseudogenes <- unique(genes$gene_id[str_detect(genes$gene_biotype, "pseudogene")])
genes <- genes[!genes$gene_id %in% pseudogenes]
# exclude duplicated loci
genes_granges <- genes[!paste0(seqnames(genes), ranges(genes)) %>% duplicated,]
# tidy up
rm(genes, pseudogenes)

#intergenic
# set seqlengths of genes_granges
chrom_lengths <- chrom_lengths <- read_delim("genome/chrom_lengths.txt", delim = "\t")
seqlengths(genes_granges) <- chrom_lengths$length
intergenic_granges <- setdiff(as(seqinfo(genes_granges), "GRanges"), genes_granges, ignore.strand = TRUE)
#tidy up
rm(chrom_lengths)

#promoters
promoters_granges <- unique(promoters(genes_granges))

# repeats
repeats <- read_delim("monDom5_UCSC_RepeatMasker_Xshifted_ClassFamily.txt",
                      delim = "\t",
                      col_names = TRUE)
repeats_granges <- GRanges(seqnames = repeats$chrom ,
                           IRanges(start = repeats$start,
                                   end = repeats$end),
                           strand = repeats$strand)

# cgis
cgis_granges <- readRDS("mondom5_CGI_Xshifted.RDS")

# list of all features
features <- list(repeats_granges, genes_granges, intergenic_granges, promoters_granges, cgis_granges)
features_names <- c( "repeatRegions", "geneRegions", "intergenicRegions", "promoterRegions", "cgiRegions")


# variables

cutoff <- 20


# load data
# make inputs for methRead function
# sample names
samples_emb <- as.list(c("sperm", "oocyte", 
                         "7.5_ED", "7.5_TE"))

samples_ad <- as.list(c("brain", "liver", "spleen")) 

samples <- c(samples_emb, samples_ad)

# files to import from
files_emb <- as.list(paste0(dataDir,
                            samples_emb,
                            "_destranded_pooled_df_filt5.csv"))

files_ad <- as.list(paste0(adultDir,
                           samples_ad,
                           "_destranded_pooled_df_filt5.csv"))

files <- (c(files_emb, files_ad))

# 'treatment' vector
condition <- c(1:length(files))

# read data
meth_dat <- methRead(location = files, 
                     sample.id = samples, 
                     assembly = "MonDom5", 
                     pipeline = list(fraction = TRUE, 
                                     chr.col = 1, 
                                     start.col = 2, 
                                     end.col = 3, 
                                     coverage.col = 5, 
                                     strand.col = 4, 
                                     freqC.col = 8
                     ), 
                     context = "CpG", 
                     resolution = "base", 
                     treatment = condition, 
                     mincov = 5)

# differential methylation
# unite to keep only rows present in both
# define sample comparisons
pair1 <- c(1, 3, 3, 3, 3)
pair2 <- c(2, 4, 5, 6, 7)



counts <- matrix(NA, ncol = length(pair1), nrow = 4, 
                 dimnames = list(c("hypo",
                                   "same",
                                   "hyper",
                                   "total"), 
                                 paste(pair1,
                                       pair2,
                                       sep = "->")
                 )
)
counts_catcher <- list()
diffs <- list()
diffs_catcher <- list()
for(z in 1:length(features)){
  

for(i in 1:length(pair1)) {
  a <- pair1[i]
  b <- pair2[i]
  pair_list <- new("methylRawList", 
                   meth_dat@.Data[c(a, b)], 
                   treatment = c(a, b)
  )
  pair_feat <- selectByOverlap(pair_list, features[[z]])
  pair_united <- methylKit::unite(pair_feat)
  pair_diffMeth <- calculateDiffMeth(pair_united)
  pair_prop <- diffMethPerChr(pair_diffMeth, plot = F, meth.cutoff = cutoff)
  pair_counts <- data.frame(
    hypo = pair_prop$diffMeth.all$percentage.of.hypomethylated, 
    hyper = pair_prop$diffMeth.all$percentage.of.hypermethylated
  )
  pair_counts$same <- 100 - pair_counts$hypo - pair_counts$hyper
  
  counts[, i] <- c(unlist(pair_counts)[c(1, 3, 2)], nrow(pair_united))
  diffs[[i]] <- pair_diffMeth
}

names(diffs) <- paste(pair1, pair2, sep = "->")
counts_catcher[[z]]  <- counts
diffs_catcher[[z]] <- diffs
}

names(counts_catcher) <- features_names
names(diffs_catcher) <- features_names


# save loop output to file
setwd(dataDir)
saveRDS(counts_catcher, file = paste0("counts_catcher_diffMeth_",
                                      cutoff,
                                      "_overlaps.RDS"))
saveRDS(diffs_catcher, file = paste0("diffs_catcher_diffMeths_",
                                     cutoff,
                                     "_overlaps.RDS"))


# reorganise
counts <- ldply(counts_catcher, rbind, .id = "region") %>%
  mutate(change = rep(c("hypo",
                        "same",
                        "hyper",
                        "total"), 5)) %>%
  # explicitly set factor levels to control order on plot
  mutate(region = factor(region,
                         levels = c("cgiRegions",
                                    "promoterRegions",
                                    "geneRegions",
                                    "intergenicRegions",
                                    "repeatRegions")))



# main figure
setwd(plotDir)


# labels
labels_plot = c("1->2" = "sperm_oocyte",
                "3->4" =  "E7.5_ED_TE",
                "3->5" = "E7.5_ED_brain",
                "3->6" = "E7.5_ED_liver",
                "3->7" = "E7.5_ED_spleen")


for(i in 1:length(labels_plot)){
counts %>%
  as.data.frame() %>%
  gather(-change, - region, key = comparison, value = perc) %>%
  filter(!change %in% c("same", "total")) %>%
  mutate(comparison_fact = factor(comparison,
                                  levels = c("1->2",
                                             "3->4",
                                             "3->5",
                                             "3->6",
                                             "3->7"))) %>%
    filter(comparison_fact == names(labels_plot)[i]) %>%
  ggplot(aes(x = "a",
             y = perc,
             colour = change,
             fill = change))+
    geom_col(size = 1.2,
             alpha = 0.6)+
  scale_fill_manual(values = c("#C02324", "#3A5FCD"),
                    labels = c("higher","lower"),
                    name = "methylation\ndifference")+
    scale_colour_manual(values = c("#861819", "#22397b"),
                        guide = FALSE)+
  facet_grid(.~region)+
  ylab("CpGs (%)")+
  xlab("genomic feature type")+
  theme_bw()+
  ylim(0,10)+
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18))
  setwd(plotDir)
ggsave(paste0(labels_plot[i], "_diff_meth_bar.pdf"), width = 12, height = 2.5)
}

# density plots
# per sample pair
# get methylation data for each feature type 
# store in list

pairs <- list()
pairs_catcher <- list()
for(i in 1:length(pair1)){
for(z in 1:length(features_names)){
  a <- pair1[i]
  b <- pair2[i]
  pair_list <- new("methylRawList", 
                   meth_dat@.Data[c(a, b)], 
                   treatment = c(a, b)
  )
  pair_feat <- selectByOverlap(pair_list, features[[z]])
  pair_united <- methylKit::unite(pair_feat)
  pairs[[z]] <- getData(pair_united) %>%
    mutate(ratio1 = numCs1/coverage1,
           ratio2 = numCs2/coverage2)
  
}
  names(pairs) <- features_names
  pairs_catcher[[i]] <- ldply(pairs, rbind, .id = "region") %>% data_frame()
}
names(pairs_catcher) <- paste0(pair1, "->", pair2)

# reformat list produced by loop to dataframe
compares <- ldply(pairs_catcher,
                  rbind,
                  .id = "comparison") %>%
  data_frame() %>%
  pivot_longer(starts_with("ratio"),
               names_to = "sample",
               values_to = "ratio") %>%
  dplyr::select(c("region",
                  "comparison",
                  "chr",
                  "start",
                  "ratio",
                  "sample")) %>%
  # explicitly set factor levels to control order on plot
  mutate(region = factor(region, levels = c("cgiRegions",
                                            "promoterRegions",
                                            "geneRegions",
                                            "intergenicRegions",
                                            "repeatRegions")))


# colours
gamete_colours <- c("#9467BD",
                  "#BCBD22")
E7.5_colours <- c("#f66e6e",
                  "#8B2323")
brain_colours <- c("#f66e6e",
                   "#27408B")
liver_colours <- c("#f66e6e",
                   "#3A5FCD")
spleen_colours <- c("#f66e6e",
                    "#436EEE")


colours <- list(gamete_colours,
                E7.5_colours,
                brain_colours,
                liver_colours,
                spleen_colours)

# make density plot per sample pair
for(i in 1:length(labels_plot)){
  compares %>% 
    filter(comparison == names(labels_plot)[i]) %>%
    ggplot(aes(x = ratio,
               y = 1,
               fill = sample,
               colour = sample,
               group = sample))+
    geom_density_ridges(alpha = 0.2,
                        size = 1.2,
                        show.legend = T)+
    facet_grid(.~region)+
    scale_fill_manual(values = colours[[i]], )+
    scale_colour_manual(values = colours[[i]])+
    theme_bw()+
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18))+
    xlim(-0.2, 1.2)
  setwd(plotDir)
  ggsave(paste0(labels_plot[i], "fig3_densities.pdf"), width = 12, height = 2)
}

  


  
