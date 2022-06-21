# BJL
# BS-seq data
#methylation at region of interest 

# libraries
library(tidyverse)
library(methylKit)
library(GenomicFeatures)
library(grid)
library(rtracklayer)


# directories
dataDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/methylKit_out/sex_specific"
annotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/annotations/"
plotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/plots/"


# define genomic region of interest
setwd(annotDir)
# regions of interest
regions <- read.delim("rsx_regions_mondom5_shifted.txt")
# make GRanges object
rois <- GRanges(seqnames = regions$chr, ranges = IRanges(regions$start, regions$end), strand = regions$strand, name = regions$name)
# defines region of interest plotted in this script
x = 11


# gene track
gtf <- import(file.path(annotDir, "mondom5.97_Xshifted_Y_DNMT1.gtf")) 
# extract gene annotations for region of interest
genes <- subsetByOverlaps(gtf, rois[x])  
# keep only exons
genes <- genes[genes$type == "exon",] 

# manually coded dependent on number of genes in roi
genes <- as.data.frame(genes)

# fudge a plot trick variable for positioning geom_rect for each gene
counter_min = 0.05
counter_max = 0.15

for(i in 1:length(unique(genes$gene_id))){
  genes$plot_trick_min[i] = counter_min
  genes$plot_trick_max[i] = counter_max
  counter_min = counter_max + 0.03
  counter_max = counter_max + 0.13
}


# loop through stages
# import meth data
# calculate meth ratio
# plot male vs female

stages <- c("1.5_embryo", "2.5_embryo", "3.5_embryo",
            "4.5_embryo", "5.5_embryo", "6.5_embryo",
            "7.5_embryo")

for(t in 1:length(stages)){
  
  stage <- stages[t]
  

# import objects for methRead
samples <- as.list(paste0(stage, c("_female", "_male")))
files <- as.list(paste0(samples, "_destranded_pooled_df_filt5.csv"))
condition <- c(1:length(samples))

# import data
setwd(dataDir)
meth_data <- methRead(location = files, 
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
                      mincov = 5
)



# make a list containing a methylation dataframe for sample pair in samples vector
meth_store <- list()
for(i in 1:length(samples)){
  meth <- as(meth_data[[i]], "GRanges")
  roi_meth <- subsetByOverlaps(meth, rois[x])
  meth_store[[i]] <- as.data.frame(roi_meth) %>%
    dplyr::select(seqnames, start, end, numCs, coverage) %>%
    mutate(ratio = numCs/coverage)
}
names(meth_store) <- samples

# bind to one dataframe
meth_df <- dplyr::bind_rows(meth_store, .id = 'source') %>%
  separate(source, into = c("stage", "sex"), sep = "_embryo_")

plot_cols_1 <- c("#AEC7E8",
                "#98DF8A",
                "#DBDB8D",
                "#FFBB78",
                "#FF9896",
                "#E377C2",
                "#CD3333")
plot_cols_2 <- c("#2c5d9e",
                 "#378d26",
                 "#88882b",
                 "#bb5d00",
                 "#ca0300",
                 "#8f1d6c",
                 "#661919")
                 
  spots_plot <- ggplot()+
    geom_point(data = meth_df,
               aes(x = start,
                   y = ratio,
               colour = sex),
               size = 2)+
    theme_classic()+
    ylab("CpG methylation")+
    # set x axis as chromosome name
    xlab(paste0(unique(regions$chr), " (Mb)"))+
    # set y scale with room for drawing genes under plot
    scale_y_continuous(breaks = c(0, 0.5, 1),
                       labels = c("0", "0.5", "1"),
                       limits = c(NA, 1))+
    # set x scale as position on chromosome in Mb
    scale_x_continuous(breaks = c(regions$start[x]-100,
                                  (regions$start[x] + regions$end[x])/2,
                                  regions$end[x]+100),
                       labels = c((regions$start[x]-100)/1000000,
                                  ((regions$start[x] + regions$end[x])/2)/1000000,
                                  (regions$end[x]+100)/1000000),
                       limits = c(regions$start[x]-100,
                                  regions$end[x]+100),
                       expand = c(0, 0))+
    scale_colour_manual(values = c(plot_cols_1[t], plot_cols_2[t]))+
    # make an hline at zero to account for the negative ylim set above
    #geom_hline(yintercept = 0)+
    # draw gene positions under plot
    # intron lines
    geom_segment(data = genes,
                 aes(y = -(plot_trick_min + plot_trick_max)/2,
                     yend = -(plot_trick_min + plot_trick_max)/2),
                 xend = regions$end[x]+100,
                 x = regions$start[x]-100)+
    # exon boxes
    geom_rect(data = genes,
              aes(xmin = start,
                  xmax = end,
                  ymin = -plot_trick_min,
                  ymax = -plot_trick_max),
              fill = "black") +
    theme(axis.ticks = element_blank(),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          plot.margin=unit(c(1,1.5,1,1),"cm"),
          legend.position = "top")
  setwd(plotDir)
  spots_plot
  ggsave(paste0(stage, "_", regions$name[x], "_spotplot.pdf"), width = 7, height = 3, useDingbats = FALSE)

}
