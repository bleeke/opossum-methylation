# BJL
# histograms of mouse control BSseq data

# libraries
library(tidyverse)
library(methylKit)
library(scales)


# directories
dataDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/mouse_bsseq/methylKit_out"
plotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/plots"

# variables and function

# y breaks function
breaks_fun = function(x) {
  c(0, max(x)*0.25, max(x)*0.5, max(x)*0.75, max(x))
}

# top labels
my_labels = function(x){
  c(rep("", times = length(x)-1), formatC((x[length(x)]/1000000), format = "f", digits = 1))
}


# get data: 'bulk' 
# create list of samples
samples <- as.list(c("sperm",
                     "E3.5_embryo",
                     "brain"))
condition <- c(1, 2, 3)
setwd(dataDir)
files <- as.list(paste0(samples,
                        "_destranded_pooled_df_filt5.csv"))
data_bulk <- methRead(location = files,
                      sample.id = samples,
                      assembly = "mm10",
                      pipeline = list(fraction = TRUE,
                                      chr.col = 1,
                                      start.col = 2,
                                      end.col = 3,
                                      coverage.col = 5,
                                      strand.col = 4,
                                      freqC.col = 8),
                      context = "CpG",
                      resolution = "base",
                      treatment = condition,
                      mincov = 5)
# apply unite function to make methylBase format but keep all sites
data_bulk <- unite(data_bulk, min.per.group = 0L)
# calculate percentage methylation
pc_data_bulk <- percMethylation(data_bulk) %>%
  data.frame() %>%
  rownames_to_column(var = "position")

# tidy up
rm(data_bulk)

# plot

# colours
cols <- c("#F0E442",
          "#CC79A7",
          "#000000")



# histogram

hist_plot <- pc_data_bulk %>%
    dplyr::rename(E3.5 = E3.5_embryo) %>%
    pivot_longer(!position,
                 names_to = "sample",
                 values_to = "pcMeth") %>%
    mutate(sample = factor(sample,
                           levels = c("sperm",
                                      "E3.5",
                                      "brain"))) %>%
    ggplot(aes(x = pcMeth,
               fill = sample))+
    # specify breaks for histogram bins
    geom_histogram(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100))+
    theme_classic()+
    xlab("CpG methylation (%)")+
    ylab("frequency")+
    # facet on sample and allow y axis scale to be free
    facet_wrap(.~sample,
               scales = "free_y",
               nrow = 1)+
    # set scientific notation and use breaks_fun to set y label 
    scale_y_continuous(labels = my_labels,
                       breaks = breaks_fun,
                       expand = expansion(mult = c(0.05,0.1)))+
    scale_x_continuous(breaks = c(0, 50, 100),
                       labels = c(0, 0.5, 1))+
    # set colours
    scale_fill_manual(values = cols)+
    # set text sizes, y label position inside margin, remove legend
    theme(axis.title = element_blank(),
          axis.ticks.length = unit(-0.25, "cm"),
          axis.text.x = element_text(size = 18,
                                     vjust = -2),
          plot.margin = unit(c(15, 5.5, 15, 5.5), "points"),
          axis.text.y = element_text(size = 18,
                                     margin = margin(l = 20, r = -33)),
          legend.position = "none",
          strip.text = element_blank())
  
  setwd(plotDir)
  hist_plot
  ggsave(filename = paste0("mouse_methylation_histogram_.pdf"),
         width = 10,
         height = 3)

# calculate mean methylation for in-text reference in ms
  
sperm_mean <- mean(pc_data_bulk$sperm, na.rm=T)
E3.5_mean <- mean(pc_data_bulk$E3.5_embryo, na.rm=T)
brain_mean <- mean(pc_data_bulk$brain, na.rm=T)




