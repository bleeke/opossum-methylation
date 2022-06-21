# figure s2
# histograms of CpG methylation by genomic feature


# libraries
library(tidyverse)
library(scales)
library(GenomicRanges)
library(rtracklayer)



# directories
dataDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/methylKit_out/"
adultDir  <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/adult_bsseq/methylKit_out"
plotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/plots/"
annotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/annotations/"



# variables and function

# y breaks function
breaks_fun = function(x) {
  c(0, max(x)*0.25, max(x)*0.5, max(x)*0.75, max(x))
}


# top labels
my_labels = function(x){
  c(rep("", times = length(x)-1), formatC((x[length(x)]/1000000), format = "f", digits = 2))
}



# get annotations 
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


# make plots


# sample labels
sample_labs_7.5 <- c("7.5_ED" = "ED",
                     "7.5_TE" = "TE")


# colours
cols <- c("#9467BD",
          "#BCBD22",
          "#AEC7E8",
          "#98DF8A",
          "#DBDB8D",
          "#FFBB78",
          "#FF9896",
          "#E377C2",
          "#CD3333")
cols_7.5 <- c("#8B2323",
              "#f66e6e")
cols_adult <- c("#27408B",
                "#3A5FCD", 
                "#436EEE")



# embryo data
setwd(dataDir)

embryo_data <- list.files(pattern = "_destranded_pooled_df_filt") %>% 
  set_names(str_replace(.,"_d.*", "")) %>%
  map_dfr(read_delim,
          delim = "\t",
          .id = "file") %>%
  mutate(percent = ratio * 100)

data_granges <- makeGRangesFromDataFrame(embryo_data)

# embryo histograms 
for(i in 1:length(features)){
  temp_index <- findOverlaps(data_granges, features[[i]])
  tmp_data <- embryo_data[queryHits(temp_index),]
  
  hist_plot <- tmp_data %>%
    filter(!file %in% c("7.5_ED", "7.5_TE")) %>%
    mutate(sample = factor(file,
                           levels = c("sperm", "oocyte", "1.5_embryo",
                                      "2.5_embryo", "3.5_embryo", "4.5_embryo",
                                      "5.5_embryo", "6.5_embryo", "7.5_embryo"))) %>%
    ggplot(aes(x = percent,
               colour = sample))+
    # specify breaks for histogram bins
    geom_histogram( # specify breaks for histogram bins
      breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
      # make fill invisible
      alpha = 0,
      # set line size
      size = 2)+
    theme_classic()+
    xlab("CpG methylation (%)")+
    ylab("CpGs (millions")+
    # facet on sample and allow y axis scale to be free
    facet_wrap(.~sample,
               scales = "free_y",
               nrow = 1)+
    # set scientific notation and use breaks_fun to set y label 
    scale_y_continuous(labels = my_labels,
                       breaks = breaks_fun,
                       expand = expansion(mult = c(0.05, 0.1)))+
    scale_x_continuous(breaks = c(0, 50, 100),
                       labels = c(0, 0.5, 1))+
    # set colours
    scale_colour_manual(values = cols)+
    #scale_fill_manual(values = cols)+
    # set text sizes, y label position inside margin, remove legend
    theme(axis.title = element_blank(),
          axis.ticks.length = unit(-0.25, "cm"),
          axis.text.x = element_text(size = 22,
                                     vjust = -1.5),
          plot.margin = unit(c(15, 5.5, 15, 5.5), "points"),
          axis.text.y = element_text(size = 22,
                                     margin = margin(l = 20, r = -45)),
          legend.position = "none",
          strip.background = element_blank(),
          strip.text.x = element_blank())
  
  setwd(plotDir)
  hist_plot
  ggsave(filename = paste0("methylation_histogram", features_names[i], "_embryo.pdf"), width = 27, height = 3.2)
}
  
# E7.5 ED and E7.5 TE overlaid histograms  

for(i in 1:length(features)){
  temp_index <- findOverlaps(data_granges, features[[i]])
  tmp_data <- embryo_data[queryHits(temp_index),]
  
  hist_plot <- tmp_data %>%
    filter(file %in% c("7.5_ED", "7.5_TE")) %>%
    mutate(sample = factor(file,
                           levels = c("7.5_TE", "7.5_ED"))) %>%
    ggplot(aes(x = percent,
               colour = sample,
               group = sample))+
    
    geom_histogram(position = "identity",
                   # specify breaks for histogram bins
                   breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
                   # make fill invisible
                   alpha = 0,
                   # set line size
                   size = 1.2)+
    theme_classic()+
    xlab("CpG methylation (%)")+
    ylab("CpGs (millions")+
    # set scientific notation and use breaks_fun to set y label 
    scale_y_continuous(labels = my_labels,
                       breaks = breaks_fun,
                       expand = expansion(mult = c(0.05, 0.1)))+
    scale_x_continuous(breaks = c(0, 50, 100),
                       labels = c(0, 0.5, 1))+
    # set colours
    scale_colour_manual(values = cols_7.5,
                        labels = sample_labs_7.5,
                        name = "")+
    # set text sizes, y label position inside margin, remove legend
    theme(axis.title = element_blank(),
          axis.ticks.length = unit(-0.25, "cm"),
          axis.text.x = element_text(size = 22,
                                     vjust = -1.5),
          plot.margin = unit(c(15, 5.5, 15, 5.5), "points"),
          axis.text.y = element_text(size = 22,
                                     margin = margin(l = 20, r = -45)),
          legend.position = "top",
          legend.text = element_text(size = 18))
  
  setwd(plotDir)
  hist_plot
  ggsave(filename = paste0("methylation_histogram_7.5_", features_names[i], "_embryo.pdf"),  width = 3, height = 3.5)
}




# adult data
setwd(adultDir)

adult_data <- list.files(pattern = "_destranded_pooled_df_filt") %>% 
  set_names(str_replace(.,"_d.*", "")) %>%
  map_dfr(read_delim,
          delim = "\t",
          .id = "file") %>%
  mutate(percent = ratio * 100) 

data_adult_granges <- makeGRangesFromDataFrame(adult_data)


# adult histograms
for(i in 1:length(features)){
  temp_index <- findOverlaps(data_adult_granges, features[[i]])
  tmp_data <- adult_data[queryHits(temp_index),]

hist_plot <- tmp_data %>%
  mutate(sample = factor(file,
                         levels = c("brain",
                                    "liver",
                                    "spleen"))) %>%
  ggplot(aes(x = percent,
             colour = sample))+
  # specify breaks for histogram bins
  geom_histogram( # specify breaks for histogram bins
    breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
    # make fill invisible
    alpha = 0,
    # set line size
    size = 2)+
  theme_classic()+
  xlab("CpG methylation (%)")+
  ylab("CpGs (millions")+
  # facet on sample and allow y axis scale to be free
  facet_wrap(.~sample,
             scales = "free_y",
             nrow = 1)+
  # set scientific notation and use breaks_fun to set y label 
  scale_y_continuous(labels = my_labels,
                     breaks = breaks_fun,
                     expand = expansion(mult = c(0.05, 0.1)))+
  scale_x_continuous(breaks = c(0, 50, 100),
                     labels = c(0, 0.5, 1))+
  # set colours
  scale_colour_manual(values = cols_adult)+
  #scale_fill_manual(values = cols)+
  # set text sizes, y label position inside margin, remove legend
  theme(axis.title = element_blank(),
        axis.ticks.length = unit(-0.25, "cm"),
        axis.text.x = element_text(size = 22,
                                   vjust = -1.5),
        plot.margin = unit(c(15, 5.5, 15, 5.5), "points"),
        axis.text.y = element_text(size = 22,
                                   margin = margin(l = 20, r = -45)),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank())

setwd(plotDir)
hist_plot
ggsave(filename = paste0("methylation_histogram_adult_", features_names[i], ".pdf"), width = 8, height = 3.6)
}
