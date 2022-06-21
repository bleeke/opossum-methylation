# BJL
# fig1 opossum methylation histograms

# libraries
library(tidyverse)
library(scales)



# directories
dataDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/methylKit_out/"
adultDir  <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/data/adult_bsseq/methylKit_out"
plotDir <- "/camp/lab/turnerj/working/Bryony/manuscript/analysis/plots/"


# variables and function

# y breaks function
breaks_fun = function(x) {
  c(0, max(x)*0.25, max(x)*0.5, max(x)*0.75, max(x))
}


# top labels
my_labels = function(x){
  c(rep("", times = length(x)-1), formatC((x[length(x)]/1000000), format = "f", digits = 1))
}

# get data


# embryo data
setwd(dataDir)

embryo_data <- list.files(pattern = "_destranded_pooled_df_filt") %>% 
  set_names(str_replace(.,"_d.*", "")) %>%
  map_dfr(read_delim,
          delim = "\t",
          .id = "file") %>%
  mutate(percent = ratio * 100)


# adult data
setwd(adultDir)

adult_data <- list.files(pattern = "_destranded_pooled_df_filt") %>% 
  set_names(str_replace(.,"_d.*", "")) %>%
  map_dfr(read_delim,
          delim = "\t",
          .id = "file") %>%
  mutate(percent = ratio * 100) 


# make histogram

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

# plot embryo data
hist_plot <- embryo_data %>%
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
                                   margin = margin(l = 20, r = -33)),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank())

setwd(plotDir)
hist_plot
ggsave(filename = paste0("methylation_histogram_all_sites_embryo.pdf"), width = 27, height = 3.2)


# overlay of E7.5 lineages
hist_plot <- embryo_data %>%
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
                                   margin = margin(l = 20, r = -33)),
        legend.position = "top",
        legend.text = element_text(size = 18))

setwd(plotDir)
hist_plot
ggsave(filename = paste0("methylation_histogram_all_sites_lineages7.5.pdf"), width = 3, height = 3.5)


# brain - main fig
hist_plot <- adult_data %>%
  mutate(sample = factor(file)) %>%
  filter(sample == "brain") %>%
  ggplot(aes(x = percent,
             colour = sample))+
  # specify breaks for histogram bins
  geom_histogram(position = "identity",
                breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
                 alpha = 0,
                 size = 2)+
  theme_classic()+
  xlab("CpG methylation (%)")+
  ylab("CpGs (millions)")+
  # set scientific notation and use breaks_fun to set y label 
  scale_y_continuous(labels = my_labels,
                     breaks = breaks_fun,
                     expand = expansion(mult = c(0.05, 0.1)))+
  scale_x_continuous(breaks = c(0, 50, 100),
                     labels = c(0, 0.5, 1))+
  # set colours
  scale_colour_manual(values = cols_adult[1],
                      guide = "none")+
  # set text sizes, y label position inside margin, remove legend
  theme(axis.title = element_blank(),
        axis.ticks.length = unit(-0.25, "cm"),
        axis.text.x = element_text(size = 22,
                                   vjust = -1.5),
        plot.margin = unit(c(15, 5.5, 15, 5.5), "points"),
        axis.text.y = element_text(size = 22,
                                   margin = margin(l = 20, r = -33)),
        legend.position = "top",
        legend.text = element_text(size = 18))

setwd(plotDir)
hist_plot
ggsave(filename = paste0("methylation_histogram_all_sites_adultbrain.pdf"), width = 3, height = 3.2)

